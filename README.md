# MGM
Source code for mixed graphical model estimation

The method is proposed in the paper
High dimensional mixed graphical models
by
Jie Cheng, Tianxi Li, Elizaveta Levina, Ji Zhu.

The paper can be found on http://arxiv.org/abs/1304.2810.


Please send your questions and feedback to tianxili@umich.edu.



The folder MethodCode contains the core code for the method.

sepreg_weight_flex.m is the main function to fit the model.
Example.m is a simple demo of how to use it, with additional functions of data generation, performance evaluation etc.


Note that the current implementation relies on the library glmnet for lasso solution. The library can be found on http://web.stanford.edu/~hastie/glmnet_matlab. 

The folder PaperExamples contains the scripts for simulation and real data examples in the paper.

The comparison with overlapping group lasso needs the implementation
of SLEP_package_4.1, which can be found on
http://www.yelab.net/software/SLEP/

The implementation of Lee & Hastie can be found on Jason Lee’s
website. We simply add a wrapper for it. So to run that method, one has to make sure all of the source files of their implementation are available and in the searching path. Note that the library TFOCS is needed for it. Please check http://www-bcf.usc.edu/~lee715/learningmgm.html for details.


The directories and paths may have to be changed accordingly.
