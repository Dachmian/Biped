using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Symbolics;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
using System.IO;
using System.Globalization;
using NagLibrary;
using MathNet.Numerics.Integration;
using NLoptNet;

namespace BipedRobot
{
    class SLSQPCoM
    {
        private BRgait _gait;
        private BRParameters _param;
        double _P0;
        double _P6;
        double _P12;

        private double[] _currentParameterValues;
        private static double _step = 1e-6;

        public delegate double func(Expression exp);

        public SLSQPCoM(Biped biped, int numOfParams)
        {
            _gait = biped.gait;
            _param = biped.param;
            BRVHC vhc = biped.gait.vhc;
            _currentParameterValues = new double[3 * numOfParams];
            while (true) {
                Random rndm = new Random();
                //map from 0 to 1 to 0.07-0.157
                _P0 = (((rndm.NextDouble() - 0.0) * (0.157 - 0.07)) / (1 - 0) + 0.07);
                _P6 = (((rndm.NextDouble() - 0.0) * (0.157 - 0.07)) / (1 - 0) + 0.07);
                _P12 = -_P0;
                double[] theta = thetaRange(_P0, _P6, _P12);
                double[] x = new double[3 * numOfParams - 3];
                for(int i = 0; i < 3 * numOfParams - 3; i++)
                {
                    x[i] = (((rndm.NextDouble() - 0.0) * (5 - (-5))) / (1 - 0) + (-5));
                }
                double[,] c = new double[,] { { theta[1], theta[1], theta[1], theta[1], theta[1], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2 * _P0 }, { 0, 0, 0, 0, 0, theta[1], theta[1], theta[1], theta[1], theta[1], 0, 0, 0, 0, 0, 0 }, { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, theta[1], theta[1], theta[1], theta[1], theta[1], -2 * _P12 } };
                int[] ct = new int[] { 1, 1, 1};
                alglib.minbleicstate state;
                alglib.minbleicreport rep;

                //
                // These variables define stopping conditions for the optimizer.
                //
                // We use very simple condition - |g|<=epsg
                //
                double epsg = 0.001;
                double epsf = 0;
                double epsx = 0;
                int maxits = 0;
                double diffstep = 1.0e-6;
                //
                // Now we are ready to actually optimize something:
                // * first we create optimizer
                // * we add boundary constraints
                // * we tune stopping conditions
                // * and, finally, optimize and obtain results...
                //
                alglib.minbleiccreatef(x, diffstep, out state);
                alglib.minbleicsetlc(state, c, ct);
                alglib.minbleicsetcond(state, epsg, epsf, epsx, maxits);
                alglib.minbleicoptimize(state, BLEICobjFunction, null, null);
                alglib.minbleicresults(state, out x, out rep);

                //
                // ...and evaluate these results
                //
                System.Console.WriteLine("{0}", rep.terminationtype); // EXPECTED: 4
                System.Console.WriteLine("{0}", alglib.ap.format(x, 2)); // EXPECTED: [-1,1]
                System.Console.ReadLine();
                double test1 = _P0 - (_P12 + theta[1] * x[10] + theta[1] * x[11] + theta[1] * x[12] + theta[1] * x[13] + theta[1] * x[14]);
                double test2 = _P6 - (_P6 + theta[1] * x[5] + theta[1] * x[6] + theta[1] * x[7] + theta[1] * x[8] + theta[1] * x[9]);
                double test3 = _P12 - (_P0 + theta[1] * x[0] + theta[1] * x[1] + theta[1] * x[2] + theta[1] * x[3] + theta[1] * x[4]);
                if ((Math.Abs(test1) + Math.Abs(test2) + Math.Abs(test3)) < 0.001)
                {
                    _currentParameterValues = new double[] { _P0, x[0], x[1], x[2], x[3], x[4], _P6, x[5], x[6], x[7], x[8], x[9], _P12, x[10], x[11], x[12], x[13], x[14] };
                    break;
                }
            }
            
        }


        public void BLEICobjFunction(double[] x, ref double func, object obj)
        {
            double[] theta = thetaRange(_P0, _P6, _P12);
            func = Math.Pow((_P0 - (_P12 + theta[1] * x[10] + theta[1] * x[11] + theta[1] * x[12] + theta[1] * x[13] + theta[1] * x[14])),2) + Math.Pow((_P6 - (_P6 + theta[1] * x[5] + theta[1] * x[6] + theta[1] * x[7] + theta[1] * x[8] + theta[1] * x[9])), 2) + Math.Pow((_P12 - (_P0 + theta[1] * x[0] + theta[1] * x[1] + theta[1] * x[2] + theta[1] * x[3] + theta[1] * x[4])), 2);
        }



        public void run()
        {
            using (var solver = new NLoptSolver(NLoptAlgorithm.LD_SLSQP, 18, 0.0001, 100))
            {
                solver.SetLowerBounds(new[] { -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0 });
                solver.SetUpperBounds(new[] { 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 });

                solver.SetMinObjective(objfun);
                solver.AddLessOrEqualZeroConstraint(alphaMaxConstraint, 0.001);
                solver.AddLessOrEqualZeroConstraint(dthetaConstraint, 0.00001);
                solver.AddEqualZeroConstraint(impactConstraint1, 0.001);
                solver.AddEqualZeroConstraint(impactConstraint2, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint1, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint2, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint3, 0.001);

                double? finalScore;

                var initialValue = _currentParameterValues;
                var result = solver.Optimize(initialValue, out finalScore);

            }

        }
        public void runAlphaAndDtheta()
        {
            using (var solver = new NLoptSolver(NLoptAlgorithm.LD_SLSQP, 18, 0.0001, 100))
            {
                solver.SetLowerBounds(new[] { -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0 });
                solver.SetUpperBounds(new[] { 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 });

                solver.SetMinObjective(objfunAlphaAndDtheta);
                solver.AddEqualZeroConstraint(geometryConstraint1, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint2, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint3, 0.001);

                double? finalScore;

                var initialValue = _currentParameterValues;
                var result = solver.Optimize(initialValue, out finalScore);

            }

        }
        public void runImpact()
        {
            using (var solver = new NLoptSolver(NLoptAlgorithm.LD_SLSQP, 18, 0.0001, 100))
            {
                solver.SetLowerBounds(new[] { -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0, -100.0 });
                solver.SetUpperBounds(new[] { 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0 });

                solver.SetMinObjective(objfunImpact);
                solver.AddLessOrEqualZeroConstraint(alphaMaxConstraint, 0.001);
                solver.AddLessOrEqualZeroConstraint(dthetaConstraint, 0.00001);
                solver.AddEqualZeroConstraint(geometryConstraint1, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint2, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint3, 0.001);

                double? finalScore;

                var initialValue = _currentParameterValues;
                var result = solver.Optimize(initialValue, out finalScore);

            }

        }

        public double objfun(double[] p, double[] grad)
        {
            int len = p.Length / 3;
            double[] theta = thetaRange(p[0], p[6], p[12]);
            double thetaMin = theta[0];
            double thetaMax = theta[1];

            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

            }
            double[] values = evalTorques(thetaMin, thetaMax);

            if (grad != null)
            {
                double[] gradValues;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValues = evalTorques(thetaMin, thetaMax);
                    grad[i] = (Math.Pow(gradValues[0], 2) + Math.Pow(gradValues[1], 2) - (Math.Pow(values[0], 2) + Math.Pow(values[1], 2))) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValues = evalTorques(thetaMin, thetaMax);
                    grad[i] = (Math.Pow(gradValues[0], 2) + Math.Pow(gradValues[1], 2) - (Math.Pow(values[0], 2) + Math.Pow(values[1], 2))) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValues = evalTorques(thetaMin, thetaMax);
                    grad[i] = (Math.Pow(gradValues[0], 2) + Math.Pow(gradValues[1], 2) - (Math.Pow(values[0], 2) + Math.Pow(values[1], 2))) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return Math.Pow(values[0], 2) + Math.Pow(values[1], 2);

        }
        public double objfunAlphaAndDtheta(double[] p, double[] grad)
        {
            int len = p.Length / 3;
            double[] theta = thetaRange(p[0], p[6], p[12]);
            double thetaMin = theta[0];
            double thetaMax = theta[1];
            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

            }
            double alpha = evaluateAlphaConstraint(thetaMin, thetaMax);
            double dthetaMin = evaluateDthetaConstraint(thetaMin, thetaMax);
            if (grad != null)
            {
                double alphaVal;
                double dthetaMinVal;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    alphaVal = evaluateAlphaConstraint(thetaMin, thetaMax);
                    dthetaMinVal = evaluateDthetaConstraint(thetaMin, thetaMax);
                    grad[i] = (alphaVal + (-1) * dthetaMinVal - (alpha + (-1) * dthetaMin)) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    alphaVal = evaluateAlphaConstraint(thetaMin, thetaMax);
                    dthetaMinVal = evaluateDthetaConstraint(thetaMin, thetaMax);
                    grad[i] = (alphaVal + (-1) * dthetaMinVal - (alpha + (-1) * dthetaMin)) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    alphaVal = evaluateAlphaConstraint(thetaMin, thetaMax);
                    dthetaMinVal = evaluateDthetaConstraint(thetaMin, thetaMax);
                    grad[i] = (alphaVal + (-1) * dthetaMinVal - (alpha + (-1) * dthetaMin)) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return alpha + (-1) * dthetaMin;

        }
        public double objfunImpact(double[] p, double[] grad)
        {
            int len = p.Length / 3;
            double[] theta = thetaRange(p[0], p[6], p[12]);
            double thetaMin = theta[0];
            double thetaMax = theta[1];

            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

            }
            double impactValue1 = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactSecondLine(thetaMin, thetaMax);
            double impactValue2 = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactThirdLine(thetaMin, thetaMax);
            if (grad != null)
            {
                double gradValue1;
                double gradValue2;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue1 = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactSecondLine(thetaMin, thetaMax);
                    gradValue2 = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactThirdLine(thetaMin, thetaMax);
                    grad[i] = (Math.Pow(gradValue1, 2) + Math.Pow(gradValue2, 2) - (Math.Pow(impactValue1, 2) + Math.Pow(impactValue2, 2))) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue1 = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactSecondLine(thetaMin, thetaMax);
                    gradValue2 = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactThirdLine(thetaMin, thetaMax);
                    grad[i] = (Math.Pow(gradValue1, 2) + Math.Pow(gradValue2, 2) - (Math.Pow(impactValue1, 2) + Math.Pow(impactValue2, 2))) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue1 = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactSecondLine(thetaMin, thetaMax);
                    gradValue2 = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactThirdLine(thetaMin, thetaMax);
                    grad[i] = (Math.Pow(gradValue1, 2) + Math.Pow(gradValue2, 2) - (Math.Pow(impactValue1, 2) + Math.Pow(impactValue2, 2))) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return Math.Pow(impactValue1, 2) + Math.Pow(impactValue2, 2);

        }


        public double geometryConstraint1(double[] p, double[] grad)
        {
            double[] theta = thetaRange(p[0], p[6], p[12]);
            double thetaMin = theta[0];
            double thetaMax = theta[1];

            if (grad != null)
            {
                grad[0] = 1;
                grad[1] = 0;
                grad[2] = 0;
                grad[3] = 0;
                grad[4] = 0;
                grad[5] = 0;
                grad[6] = 0;
                grad[7] = 0;
                grad[8] = 0;
                grad[9] = 0;
                grad[10] = 0;
                grad[11] = 0;
                grad[12] = -thetaMax;
                grad[13] = -thetaMax;
                grad[14] = -thetaMax;
                grad[15] = -thetaMax;
                grad[16] = -thetaMax;
                grad[17] = -thetaMax;
            }
            return p[0] - (p[12] + p[13] + p[14] + p[15] + p[16] + p[17]);
        }
        public double geometryConstraint2(double[] p, double[] grad)
        {
            double[] theta = thetaRange(p[0], p[6], p[12]);
            double thetaMin = theta[0];
            double thetaMax = theta[1];
            if (grad != null)
            {
                grad[0] = 0;
                grad[1] = 0;
                grad[2] = 0;
                grad[3] = 0;
                grad[4] = 0;
                grad[5] = 0;
                grad[6] = 0;
                grad[7] = -thetaMax;
                grad[8] = -thetaMax;
                grad[9] = -thetaMax;
                grad[10] = -thetaMax;
                grad[11] = -thetaMax;
                grad[12] = 0;
                grad[13] = 0;
                grad[14] = 0;
                grad[15] = 0;
                grad[16] = 0;
                grad[17] = 0;
            }
            return p[6] - (p[6] + p[7] + p[8] + p[9] + p[10] + p[11]);
        }
        public double geometryConstraint3(double[] p, double[] grad)
        {
            double[] theta = thetaRange(p[0], p[6], p[12]);
            double thetaMin = theta[0];
            double thetaMax = theta[1];
            if (grad != null)
            {
                grad[0] = -thetaMax;
                grad[1] = -thetaMax;
                grad[2] = -thetaMax;
                grad[3] = -thetaMax;
                grad[4] = -thetaMax;
                grad[5] = -thetaMax;
                grad[6] = 0;
                grad[7] = 0;
                grad[8] = 0;
                grad[9] = 0;
                grad[10] = 0;
                grad[11] = 0;
                grad[12] = 1;
                grad[13] = 0;
                grad[14] = 0;
                grad[15] = 0;
                grad[16] = 0;
                grad[17] = 0;
            }
            return p[12] - (p[0] + p[1] + p[2] + p[3] + p[4] + p[5]);
        }

        public double alphaMaxConstraint(double[] p, double[] grad)
        {
            int len = p.Length / 3;
            double[] theta = thetaRange(p[0], p[6], p[12]);
            double thetaMin = theta[0];
            double thetaMax = theta[1];

            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

            }
            double alpha = evaluateAlphaConstraint(thetaMin, thetaMax);
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateAlphaConstraint(thetaMin, thetaMax);
                    grad[i] = (gradValue - alpha) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateAlphaConstraint(thetaMin, thetaMax);
                    grad[i] = (gradValue - alpha) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateAlphaConstraint(thetaMin, thetaMax);
                    grad[i] = (gradValue - alpha) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return alpha;
        }
        public double dthetaConstraint(double[] p, double[] grad)
        {
            int len = p.Length / 3;
            double[] theta = thetaRange(p[0], p[6], p[12]);
            double thetaMin = theta[0];
            double thetaMax = theta[1];

            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

            }
            double dthetaMin = evaluateDthetaConstraint(thetaMin, thetaMax);
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateDthetaConstraint(thetaMin, thetaMax);
                    grad[i] = (gradValue - dthetaMin) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateDthetaConstraint(thetaMin, thetaMax);
                    grad[i] = (gradValue - dthetaMin) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateDthetaConstraint(thetaMin, thetaMax);
                    grad[i] = (gradValue - dthetaMin) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return (-1) * dthetaMin;
        }
        public double impactConstraint1(double[] p, double[] grad)
        {
            int len = p.Length / 3;
            double[] theta = thetaRange(p[0], p[6], p[12]);
            double thetaMin = theta[0];
            double thetaMax = theta[1];

            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

            }
            double impactValue = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactSecondLine(thetaMin, thetaMax);
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactSecondLine(thetaMin, thetaMax);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactSecondLine(thetaMin, thetaMax);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactSecondLine(thetaMin, thetaMax);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return impactValue;
        }
        public double impactConstraint2(double[] p, double[] grad)
        {
            int len = p.Length / 3;
            double[] theta = thetaRange(p[0], p[6], p[12]);
            double thetaMin = theta[0];
            double thetaMax = theta[1];

            for (int i = 0; i < len; i++)
            {
                _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = len; i < 2 * len; i++)
            {
                _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

            }
            for (int i = 2 * len; i < 3 * len; i++)
            {
                _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

            }
            double impactValue = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactThirdLine(thetaMin, thetaMax);
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactThirdLine(thetaMin, thetaMax);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactThirdLine(thetaMin, thetaMax);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(thetaMin, thetaMax) - _gait.impactThirdLine(thetaMin, thetaMax);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return impactValue;

        }
        public double[] evalTorques(double thetaMin, double thetaMax)
        {
            double[][] firstIntegral = TrapezoidalSum2.calculateFirstIntegral(_gait.vhc.evalTwoTimesBetaDividedByAlpha, thetaMin, thetaMax, 200);
            double[][] secondIntegral = TrapezoidalSum2.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1], thetaMin, thetaMax, 200);
            int len = firstIntegral[1].Length;
            double dthetaTSquared = (-secondIntegral[1][len - 1] * Math.Exp(firstIntegral[1][len - 1])) / (1 - Math.Exp(firstIntegral[1][len - 1]) * Math.Pow(_gait.impactSecondLine(thetaMin, thetaMax), 2));
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(thetaMin, thetaMax), 2);
            double torque1Max = 0;
            double torque2Max = 0;
            double torque1Diff = 0;
            double torque2Diff = 0;
            double[] values;
            if (double.IsNaN(dthetaTSquared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] { torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsPositiveInfinity(dthetaTSquared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] { torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsNaN(dtheta0Squared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] { torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else if (double.IsPositiveInfinity(dtheta0Squared))
            {
                torque1Max = (Math.Pow(10, 300));
                torque2Max = (Math.Pow(10, 300));
                torque1Diff = (Math.Pow(10, 300));
                torque2Diff = (Math.Pow(10, 300));
                values = new double[] { torque1Max, torque2Max, torque1Diff, torque2Diff };
            }
            else
            {
                double[,] THETA = new double[2, len];
                THETA[0, 0] = thetaMin;
                THETA[1, 0] = dtheta0Squared;

                double torque1 = 0;
                double torque2 = 0;
                double theta = thetaMin;
                double dthetaSquared = 0;
                double ddtheta = 0;
                for (int i = 0; i < len; i += 5)
                {
                    theta = secondIntegral[0][i];
                    dthetaSquared = -secondIntegral[1][i] * Math.Exp(firstIntegral[1][i]) + Math.Exp(firstIntegral[1][i]) * dtheta0Squared;
                    ddtheta = (-_gait.vhc.evalBeta(theta) * dthetaSquared - _gait.vhc.evalGamma(theta)) / _gait.vhc.evalAlpha(theta);
                    torque1 = Math.Abs(_gait.vhc.evalAlpha1(theta) * ddtheta + _gait.vhc.evalBeta1(theta) * dthetaSquared + _gait.vhc.evalGamma1(theta));
                    torque2 = Math.Abs(_gait.vhc.evalAlpha3(theta) * ddtheta + _gait.vhc.evalBeta3(theta) * dthetaSquared + _gait.vhc.evalGamma3(theta));
                    if (torque1 > torque1Max)
                    {
                        torque1Max = torque1;
                    }
                    if (torque2 > torque2Max)
                    {
                        torque2Max = torque2;
                    }
                }
                double ddtheta0 = (-_gait.vhc.evalBeta(thetaMin) * dtheta0Squared - _gait.vhc.evalGamma(thetaMin)) / _gait.vhc.evalAlpha(thetaMin);
                torque1Diff = (_gait.vhc.evalAlpha1(thetaMin) * ddtheta0 + _gait.vhc.evalBeta1(thetaMin) * dtheta0Squared + _gait.vhc.evalGamma1(thetaMin)) -
                    (_gait.vhc.evalAlpha3(thetaMax) * ddtheta + _gait.vhc.evalBeta3(thetaMax) * dthetaTSquared + _gait.vhc.evalGamma3(thetaMax));

                torque2Diff = _gait.vhc.evalAlpha3(thetaMin) * ddtheta0 + _gait.vhc.evalBeta3(thetaMin) * dtheta0Squared + _gait.vhc.evalGamma3(thetaMin) -
                    (_gait.vhc.evalAlpha1(thetaMax) * ddtheta + _gait.vhc.evalBeta1(thetaMax) * dthetaTSquared + _gait.vhc.evalGamma1(thetaMax));
                values = new double[] { torque1Max, torque2Max, torque1Diff, torque2Diff };

            }
            return values;
        }

        public double[] thetaRange(double P0, double P6, double P12)
        {
            double a = -(_param.m1 * _param.l1 + _param.m2 * _param.L1 + _param.m3 * _param.L1) / (_param.m1 + _param.m2 + _param.m3);
            double b = -(_param.m2 * _param.l2) / (_param.m1 + _param.m2 + _param.m3);
            double c = -(_param.m3 * _param.l3) / (_param.m1 + _param.m2 + _param.m3);
            double thetaMin = a * P0 + b * P6 + c * P12;
            double thetaMax = a * (-P0) + b * (P6) + c * (-P12);
            return new double[] { thetaMin, thetaMax };
        }
        public double evaluateAlphaConstraint(double thetaMin, double thetaMax)
        {


            double dx = 0.02;
            double theta = thetaMin;
            double val = double.MinValue;
            double alphaVal = 0;
            for (int i = 0; i < thetaMax / dx; i++)
            {
                alphaVal = _gait.vhc.evalAlpha(theta);
                if (alphaVal > val)
                {
                    val = alphaVal;
                }
                theta += dx;
            }
            return val;
        }
        public double evaluateDthetaConstraint(double thetaMin, double thetaMax)
        {
            double[][] firstIntegral = TrapezoidalSum2.calculateFirstIntegral(_gait.vhc.evalTwoTimesBetaDividedByAlpha, thetaMin, thetaMax, 200);
            double[][] secondIntegral = TrapezoidalSum2.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1], thetaMin, thetaMax, 200);
            int len = firstIntegral[1].Length;
            double dthetaTSquared = (-secondIntegral[1][len - 1] * Math.Exp(firstIntegral[1][len - 1])) / (1 - Math.Exp(firstIntegral[1][len - 1]) * Math.Pow(_gait.impactSecondLine(thetaMin, thetaMax), 2));
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(thetaMin, thetaMax), 2);
            double dthetaMin = 0;
            if (double.IsNaN(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dthetaTSquared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsNaN(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else if (double.IsPositiveInfinity(dtheta0Squared))
            {
                dthetaMin = (-Math.Pow(10, 300));
            }
            else
            {
                double[,] THETA = new double[2, len];
                THETA[0, 0] = thetaMin;
                THETA[1, 0] = dtheta0Squared;
                dthetaMin = dtheta0Squared;
                for (int i = 1; i < len; i++)
                {
                    THETA[0, i] = secondIntegral[0][i];
                    THETA[1, i] = -secondIntegral[1][i] * Math.Exp(firstIntegral[1][i]) + Math.Exp(firstIntegral[1][i]) * dtheta0Squared;
                    if (THETA[1, i] < dthetaMin)
                    {
                        dthetaMin = THETA[1, i];
                    }
                }
            }
            return dthetaMin;
        }
    }
}
