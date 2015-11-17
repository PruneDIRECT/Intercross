----------------------------------------------------------------------
PruneDIRECT - R Version

copyright (c) 2015, Behrang Mahjani, Salman Toor, Carl Nettelblad, Sverker Holmgren

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------
Here, the implementation is based on the general DIRECT implementation in the R package bootfs [1]. The DIRECT code
in that package is in turn a re-implementation in R of a Matlab code developed in [2]. We also use a function from R/QTL for
calculating the probabilites [3]. The "hyper" data is from R/QTL package [3]. 


[1] C. Bender, bootfs: Use multiple feature selection algorithms
to derive robust feature sets for two or multiclass classification
problems. R package version 1.4.2/r25, 2013, http://R-Forge.
R-project.org/projects/bootfs/

[2] D. E. Finkel, Global Optimization with the DIRECT Algorithm.,
PhD theis, North Carolina State University, Feb 2005.

[3] K.W. Broman, H. Wu, S. Sen, G.A. Churchill, R/qtl: QTL
mapping in experimental crosses. Bioinformatics, 19, pp. 889-890,
2003.

----------------------------------------------------------------------

#PruneDIRECT - V7 for intercross
## Instructions: 

Where to change: 
*1. Main file, the path for the file and loading the data, 
*2. Main file, the choromosome numbers
*3. Change the objective function in the file objectiveFB.R, the name of the function is ObjectiveF3dInt 
