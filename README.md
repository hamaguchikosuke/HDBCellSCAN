# HDBCellSCAN
Hierarchical density-based clustering method to detect ROIs in Ca imaging data
<img src="https://github.com/hamaguchikosuke/HDBCellSCAN/blob/master/CaGui/figures/HDBCellSCAN_ROIs.png" width=400px>

# I. History
HDBCellSCAN was initially developped as a customized version of Suite2P. Later, new ROI detection algorithm HDBCellSCAN was implemented by Kosuke Hamaguchi. Spike deconvolution is based on Fast-Oopsi. 
Image registration based on phase correlation is based on Suite2P code.

# II. Installation. 
**Requirement**

To run HDBCellSCAN, you need

**MATLAB 2017b** or higher

**Python 3.5 or higher** and hdbscan package 

### Install Python and hdbscan ###

1. Download Anaconda Navigator and install. 

2. Open Anaconda Navigator. Make a new environment. 
   Put a name like "CaImaging" for that environment.

3. Open terminal by right-clicking the triangle button in "CaImaging" environment.  
   
  Type 
\>> conda install -c conda-forge hdbscan

During the installation, you can find the environment location, such as 
C:\Users\hammer\AppData\Local\conda\conda\envs\CaImaging

The HDBSCAN package is in 
C:\Users\hammer\AppData\Local\conda\conda\envs\CaImaging\Lib\site-packages

Please write down this path and add python path (see How2Use_HDBCellSCAN.m) so that MATLAB can reach to HDBSCAN packages. 

4. Go to Github website 

5. Download HDBCellScan package and unzip the downloaded folder.
Hereafter, I assume that HDBCellSCAN is extracted in \<RootDir\>=C:\home\GitHub\HDBCellSCAN

6. To confirm HDBSCAN package itself is working in Python, install and launch Spyder from Anaconda. 
Open C:\home\GitHub\HDBCellSCAN\hdbscan_test20171017.py and run this code.
If you face any errors of missing package, install it from Anaconda.
 
### Run HDBCellSCAN ###
7. Start MATLAB, move to \<RootDir\>. 

8. In the command prompt of MATLAB, type
\>>pyversion

to confirm you have access to python from MATLAB. If it returns empty, you need to install Python or need to set PATH variable (see trouble shooting). 

9. Open an example file, <RootDir>\How2Use_HDBCellScan.m
 You can see how to operate HDBCellScan in this code.


# III. Trouble shooting (in Windows)

### 1. pyversion returns empty 

Recent Anaconda does not recommend to set PATH during the installation. 
You need to set PATH variable by going to 
My Computer -> Properties -> Advanced -> Environment Variables and edit "Path" variable to add Python Path i.e. 'C:\<python_library>'

Another possibility is the version compatibility.
(please check https://jp.mathworks.com/matlabcentral/answers/229501-matlab-do-not-recongnize-pyversion)

### 2. svmtrain returns an error 
Note that I use libsvm (https://www.csie.ntu.edu.tw/~cjlin/libsvm/) to run SVM instead of matlab's toolbox. 
This is in part due to historical reason. In my GitHub commit, svmtrain and svmpredict are already compiled for Windows 64bit, but not for the other environments (since I do not access to Mac/Linux environment).
Please compile them using make.m in libsvm folder.


