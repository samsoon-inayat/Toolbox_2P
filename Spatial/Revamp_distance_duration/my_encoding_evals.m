function out = my_encoding_evals(FR, stimulus, nbinsMI, nshuffles)

nbinFRSH = nbinsMI;

MI = MutualInformation(stimulus',FR');
