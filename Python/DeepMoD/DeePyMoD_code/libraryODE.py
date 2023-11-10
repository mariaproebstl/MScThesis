""" Contains the library classes that store the parameters u_t, theta"""
import numpy as np
import torch
from torch.autograd import grad
from .deepmod import Library
from typing import Tuple
from ..utils.types import TensorList


class LibraryODE(Library):
    def __init__(self, int_order: int, intercept: bool) -> None:
        """Calculates a library/feature matrix consisting of
        1) intercept (if intercept = True),
        2) growth rate,
        3) 2nd order interactions, i.e. [x1*x1, x1*x2, ..., x1*xn] and
        4) 3rd order interactions, i.e. [x1 * (x1*x1, x1*x2, ..., x1*xn, ..., xn*xn)]
        OR (in case of int_order=1)
        3) linear effects of all xi, i.e. [x1, x2, ..., xn]

        Order of terms is as in the description above

        Args:
            int_order (int): maximum order of the interactions in the library - either 2 or 3 or 1 (for only linear effects) or 4, 5 (incl. all linear, quadratic (and 3rd order) terms)
            intercept (bool): the parameter tells whether an intercept is part of the library or not.
            True means that a vector containing ones is added to the library.
        """

        super().__init__()
        self.int_order = int_order
        self.intercept = intercept

    def library(
        self, input: Tuple[torch.Tensor, torch.Tensor]
    ) -> Tuple[TensorList, TensorList]:
        """Compute the temporal derivative given by data.
            Data should have t in only column.

        Args:
            input (Tuple[torch.Tensor, torch.Tensor]): A prediction u (n_samples, n_species) and time points (n_samples, 1).

        Returns:
            [TensorList]: The thetas [(n_samples, 2 + n_outputs + ((n_outputs+1) over 2))]
            computed from the library and data.
        """
        prediction, data = input
        n_samples, n_outputs = prediction.shape
        
        theta = []
        time_deriv_list = []

        ### Construct the theta matrix

        # list all 1-dimensional combinations (1, x1, x2, ..., xn)
        comb_1D = torch.cat((torch.ones(n_samples, 1), prediction), dim = 1)
        
        # for only linear effects
        if self.int_order == 1:
            for output in np.arange(n_outputs):
                if self.intercept:
                    theta_i = comb_1D
                else: # without intercept
                    theta_i = prediction
                theta.append(theta_i)
            
        # for only 2nd order interactions
        elif self.int_order == 2:
            for output in np.arange(n_outputs):
                int_2D = torch.mul(prediction[:, output : output + 1], comb_1D)
                if self.intercept:
                    theta_i = torch.cat((torch.ones(n_samples, 1), int_2D), dim = 1)
                else: # without intercept
                    theta_i = int_2D
                theta.append(theta_i)
            
        # to include 2nd and 3rd order interactions
        elif self.int_order == 3:
            # list all 2-dimensional combinations (x1^2, x1*x2, ..., x1*xn, ..., xn^2)
            comb_2D = torch.empty(n_samples, 1)
            for i in np.arange(n_outputs):
                for j in np.arange(i, n_outputs):
                    # multiply x_i with x_j
                    interaction = torch.mul(prediction[:, i:i+1], prediction[:, j:j+1])
                    comb_2D = torch.cat((comb_2D, interaction), dim = 1)
            # remove first column of comb_2D which consists of empty values
            comb_2D = comb_2D[:, 1:]

            for output in np.arange(n_outputs):
                int_2D = torch.mul(prediction[:, output : output + 1], comb_1D)
                int_3D = torch.mul(prediction[:, output : output + 1], comb_2D)
                if self.intercept:
                    theta_i = torch.cat((torch.ones(n_samples, 1), int_2D, int_3D), dim = 1)
                else: # without intercept
                    theta_i = torch.cat((int_2D, int_3D), dim = 1)
                theta.append(theta_i)
            
        # to include linear, 2nd and 3rd order interactions
        elif self.int_order == 4:
            # list all 2-dimensional combinations (x1^2, x1*x2, ..., x1*xn, ..., xn^2)
            comb_2D = torch.empty(n_samples, 1)
            for i in np.arange(n_outputs):
                for j in np.arange(i, n_outputs):
                    # multiply x_i with x_j
                    interaction = torch.mul(prediction[:, i:i+1], prediction[:, j:j+1])
                    comb_2D = torch.cat((comb_2D, interaction), dim = 1)
            # remove first column of comb_2D which consists of empty values
            comb_2D = comb_2D[:, 1:]

            for output in np.arange(n_outputs):
                int_2D = torch.mul(prediction[:, output : output + 1], prediction)
                int_3D = torch.mul(prediction[:, output : output + 1], comb_2D)
                if self.intercept:
                    theta_i = torch.cat((comb_1D, int_2D, int_3D), dim = 1)
                else: # without intercept
                    theta_i = torch.cat((prediction, int_2D, int_3D), dim = 1)
                theta.append(theta_i)
        
        # to include all linear and 2nd-order interactions
        elif self.int_order == 5:
            for output in np.arange(n_outputs):
                int_2D = torch.mul(prediction[:, output : output + 1], prediction)
                if self.intercept:
                    theta_i = torch.cat((comb_1D, int_2D), dim = 1)
                else: # without intercept
                    theta_i = torch.cat((prediction, int_2D), dim = 1)
                theta.append(theta_i)

        else: 
            raise ValueError("int_order must be either 1, 2, 3 or 4, 5!")

        ### Construct a list of time_derivatives
        for output in np.arange(n_outputs):
            dy = grad(
                prediction[:, output],
                data,
                grad_outputs=torch.ones_like(prediction[:, output]),
                create_graph=True,
            )[0]
            time_deriv = dy[:, 0:1]
            time_deriv_list.append(time_deriv)

        return time_deriv_list, theta
    
    def get_content(
        self, input: torch.Tensor
    ) -> TensorList:
        """Give a list (chr) of which expressions are in the library.

        Args:
            input (torch.Tensor): A prediction u (n_samples, n_species).

        Returns:
            [TensorList]: The expressions of the theta [(n_samples, 2 + n_outputs + ((n_outputs+1) over 2))]
            computed from the library and data.
        """
        
        prediction = input
        n_outputs = prediction.shape[1]

        # initialize lists for string output
        x_i = [f"x{x}" for x in np.arange(1, n_outputs + 1)]
        comb_1D_chr = ["1"]
        comb_1D_chr.extend(x_i)
        
        comb_all_chr = []

        if self.int_order == 1:
            for output in np.arange(n_outputs):
                if self.intercept:
                    theta_i_chr = comb_1D_chr
                else: # without intercept
                    theta_i_chr = comb_1D_chr[1:]
                comb_all_chr.append(theta_i_chr) # all linear terms
        elif self.int_order == 2:
            for output in np.arange(n_outputs):
                int_2D_chr = [f"x{output+1}*" + x for x in comb_1D_chr]
                if self.intercept:
                    theta_i_chr = [1]
                else: # without intercept
                    theta_i_chr = []
                theta_i_chr.extend(int_2D_chr) # x_i + quadratic terms including x_i
                comb_all_chr.append(theta_i_chr)     
        elif self.int_order == 3:
            comb_2D_chr = []
            for i in np.arange(n_outputs):
                for j in np.arange(i, n_outputs):
                    comb_2D_chr.append(f"x{i+1}*x{j+1}")
            for output in np.arange(n_outputs):
                int_2D_chr = [f"x{output+1}*" + x for x in comb_1D_chr]
                int_3D_chr = [f"x{output+1}*" + x for x in comb_2D_chr]
                if self.intercept:
                    theta_i_chr = [1]
                else: # without intercept
                    theta_i_chr = []
                theta_i_chr.extend(int_2D_chr) # x_i + quadratic terms including x_i
                theta_i_chr.extend(int_3D_chr) # 3rd order terms (only including x_i)
                comb_all_chr.append(theta_i_chr)
        elif self.int_order == 4:
            comb_2D_chr = []
            for i in np.arange(n_outputs):
                for j in np.arange(i, n_outputs):
                    comb_2D_chr.append(f"x{i+1}*x{j+1}")
            for output in np.arange(n_outputs):
                int_2D_chr = [f"x{output+1}*" + x for x in comb_1D_chr[1:]]
                int_3D_chr = [f"x{output+1}*" + x for x in comb_2D_chr]
                if self.intercept:
                    theta_i_chr = [1]
                else: # without intercept
                    theta_i_chr = []
                theta_i_chr.extend(comb_1D_chr[1:]) # all linear terms
                theta_i_chr.extend(int_2D_chr) # quadratic terms (only including x_i)
                theta_i_chr.extend(int_3D_chr) # 3rd order terms (only including x_i)
                comb_all_chr.append(theta_i_chr)        
        # to include all linear and 2nd-order interactions
        elif self.int_order == 5:
            for output in np.arange(n_outputs):
                int_2D_chr = [f"x{output+1}*" + x for x in comb_1D_chr[1:]]
                if self.intercept:
                    theta_i_chr = [1]
                else: # without intercept
                    theta_i_chr = []
                theta_i_chr.extend(comb_1D_chr[1:]) # all linear terms
                theta_i_chr.extend(int_2D_chr) # quadratic terms (only including x_i)
                comb_all_chr.append(theta_i_chr)        
        else: 
            raise ValueError("int_order must be either 1, 2, 3 or 4, 5!")

        return comb_all_chr
