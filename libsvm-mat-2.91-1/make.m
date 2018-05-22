% This make.m is used under Windows

% add -largeArrayDims on 64-bit machines

mex -O -largeArrayDims -c svm.cpp
mex -O -largeArrayDims -c svm_model_matlab.c
mex -O -largeArrayDims svmtrain.c svm.obj svm_model_matlab.obj
mex -O -largeArrayDims svmpredict.c svm.obj svm_model_matlab.obj
mex -O -largeArrayDims libsvmread.c
mex -O -largeArrayDims libsvmwrite.c
