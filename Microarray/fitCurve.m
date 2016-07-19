function FitCurve = fitCurve(Dvect, Fit, f_handle,Stats,c)

switch Fit{1}

        case 'linear'
            FitCurve = c.p1*Dvect + c.p2;
        case 'exp'
            FitCurve = c.A*exp(-c.n*Dvect) + c.B;
        case 'exp_1_0'
            FitCurve = exp(-c.n*Dvect);
        case 'decay'
            FitCurve = c.A/Dvect + c.B;
            Residuals = Rvect' - FitCurve;
        case 'exp0'
            FitCurve = c.A.*exp(-c.n*Dvect);
        case 'exp1'
            FitCurve = exp(-c.n*Dvect) + c.B;
        case 'decayEta0'
            FitCurve = c.A.*Dvect.^(-c.n);
        case 'decayEta'
            FitCurve = c.A.*Dvect.^(-c.n) + c.B;
        case 'decay0'
            FitCurve = c.A./Dvect;
         
end

end