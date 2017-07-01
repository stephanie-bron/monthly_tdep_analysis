double mpfSigmaDegRescaled(double mpfSigmaDegOld, double mmueEn) {
 double x = log10(mmueEn);
 return mpfSigmaDegOld*(4.974823 -1.967809*x +0.2706778*pow(x,2.));
}
