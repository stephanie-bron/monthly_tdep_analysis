double RescaledSigma_IC86_II_III_IV_SplineMPE(double Sigma, double Energy, double Zenith) {
    // Pull correction done with a polynomial with is the result of the fit on MC (could be done with a spline fct)
    double x = log10(Energy);
    if (Zenith > 85.)  {
        return Sigma * pow(10.,(1.176 - 0.834 * pow(x,1.) + 0.241 * pow(x,2.) - 0.026 * pow(x,3.) +0.001 * pow(x,4.))); // Pb (NCh)
    }
    else {
        return Sigma * pow(10.,(1.637 - 1.265 * pow(x,1.) + 0.371 * pow(x,2.) - 0.042 * pow(x,3.) +0.002 * pow(x,4.))); // Pb (NCh)
    }  
}
