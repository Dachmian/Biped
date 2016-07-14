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
        double _q1;
        double _q2;
        double _q3;
        double _thetaMin;
        double _thetaMax;

        public double q1
        {
            get
            {
                return _q1;
            }
            set
            {
                _q1 = value;
            }
        }
        public double q2
        {
            get
            {
                return _q2;
            }
            set
            {
                _q2 = value;
            }
        }
        public double q3
        {
            get
            {
                return _q3;
            }
            set
            {
                _q3 = value;
            }
        }

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
                _q1 = (((rndm.NextDouble() - 0.0) * (0.157 - (0.07))) / (1 - 0) + (0.07));
                _q2 = -(((rndm.NextDouble() - 0.0) * (0.157 - (0.07))) / (1 - 0) + (0.07));
                _q3 = -_q1;
                thetaRange(_q1, _q2, _q3);
                double[] x = new double[3 * numOfParams];
                for(int i = 0; i < 3 * numOfParams; i++)
                {
                    x[i] = (((rndm.NextDouble() - 0.0) * (5 - (-5))) / (1 - 0) + (-5));
                }
                double[,] c = new double[,] { {1 + 1, _thetaMin + _thetaMax, Math.Pow(_thetaMin, 2) + Math.Pow(_thetaMax, 2), Math.Pow(_thetaMin, 3) + Math.Pow(_thetaMax, 3), Math.Pow(_thetaMin, 4) + Math.Pow(_thetaMax, 4), Math.Pow(_thetaMin, 5) + Math.Pow(_thetaMax, 5), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   0 }, 
                    { 0, 0, 0, 0, 0, 0, 1 - 1, _thetaMin - _thetaMax, Math.Pow(_thetaMin, 2) - Math.Pow(_thetaMax, 2), Math.Pow(_thetaMin, 3) - Math.Pow(_thetaMax, 3), Math.Pow(_thetaMin, 4) - Math.Pow(_thetaMax, 4), Math.Pow(_thetaMin, 5) - Math.Pow(_thetaMax, 5), 0, 0, 0, 0, 0, 0, 0 }, 
                    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 + 1, _thetaMin + _thetaMax, Math.Pow(_thetaMin, 2) + Math.Pow(_thetaMax, 2), Math.Pow(_thetaMin, 3) + Math.Pow(_thetaMax, 3), Math.Pow(_thetaMin, 4) + Math.Pow(_thetaMax, 4), Math.Pow(_thetaMin, 5) + Math.Pow(_thetaMax, 5),   0 },
                    {1, _thetaMin, Math.Pow(_thetaMin, 2), Math.Pow(_thetaMin, 3), Math.Pow(_thetaMin, 4), Math.Pow(_thetaMin, 5), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,   _q1 },
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, _thetaMin, Math.Pow(_thetaMin, 2), Math.Pow(_thetaMin, 3), Math.Pow(_thetaMin, 4), Math.Pow(_thetaMin, 5),   _q3 },
                    {0, 0, 0, 0, 0, 0, 1, _thetaMin, Math.Pow(_thetaMin, 2), Math.Pow(_thetaMin, 3), Math.Pow(_thetaMin, 4), Math.Pow(_thetaMin, 5), 0, 0, 0, 0, 0, 0,   _q2 },
                };
                int[] ct = new int[] { 0, 0, 0, 0, 0 ,0};
                alglib.minbleicstate state;
                alglib.minbleicreport rep;

                //
                // These variables define stopping conditions for the optimizer.
                //
                // We use very simple condition - |g|<=epsg
                //
                double epsg = 0.000001;
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
                double test1 = x[0] + _thetaMin * x[1] + Math.Pow(_thetaMin, 2) * x[2] + Math.Pow(_thetaMin, 3) * x[3] + Math.Pow(_thetaMin, 4) * x[4] + Math.Pow(_thetaMin, 5) * x[5] - (x[12] + _thetaMax * x[13] + Math.Pow(_thetaMax, 2) * x[14] + Math.Pow(_thetaMax, 3) * x[15] + Math.Pow(_thetaMax, 4) * x[16] + Math.Pow(_thetaMax, 5) * x[17]);
                double test2 = x[6] + _thetaMin * x[7] + Math.Pow(_thetaMin, 2) * x[8] + Math.Pow(_thetaMin, 3) * x[9] + Math.Pow(_thetaMin, 4) * x[10] + Math.Pow(_thetaMin, 5) * x[11] - (x[6] + _thetaMax * x[7] + Math.Pow(_thetaMax, 2) * x[8] + Math.Pow(_thetaMax, 3) * x[9] + Math.Pow(_thetaMax, 4) * x[10] + Math.Pow(_thetaMax, 5) * x[11]);
                double test3 = x[12] + _thetaMin * x[13] + Math.Pow(_thetaMin, 2) * x[14] + Math.Pow(_thetaMin, 3) * x[15] + Math.Pow(_thetaMin, 4) * x[16] + Math.Pow(_thetaMin, 5) * x[17] - (x[0] + _thetaMax * x[1] + Math.Pow(_thetaMax, 2) * x[2] + Math.Pow(_thetaMax, 3) * x[3] + Math.Pow(_thetaMax, 4) * x[4] + Math.Pow(_thetaMax, 5) * x[5]);
                if ((Math.Abs(test1) + Math.Abs(test2) + Math.Abs(test3)) < 0.01)
                {
                    _currentParameterValues = new double[] { x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], x[11], x[12], x[13], x[14], x[15], x[16], x[17] };
                    int len = _currentParameterValues.Length / 3;

                    for (int i = 0; i < len; i++)
                    {
                        _gait.vhc.phi1Parameters["P" + i.ToString()] = _currentParameterValues[i];

                    }
                    for (int i = len; i < 2 * len; i++)
                    {
                        _gait.vhc.phi2Parameters["P" + i.ToString()] = _currentParameterValues[i];

                    }
                    for (int i = 2 * len; i < 3 * len; i++)
                    {
                        _gait.vhc.phi3Parameters["P" + i.ToString()] = _currentParameterValues[i];

                    }
                    break;
                }
            }
            
        }


        public void BLEICobjFunction(double[] x, ref double func, object obj)
        {
            //func = Math.Pow(((x[0] + _thetaMin * x[1] + Math.Pow(_thetaMin, 2) * x[2] + Math.Pow(_thetaMin, 3) * x[3] + Math.Pow(_thetaMin, 4) * x[4] + Math.Pow(_thetaMin, 5) * x[5]) - (x[12] + _thetaMax * x[13] + Math.Pow(_thetaMax, 2) * x[14] + Math.Pow(_thetaMax, 3) * x[15] + Math.Pow(_thetaMax, 4) * x[16] + Math.Pow(_thetaMax, 5) * x[17])), 2) +
            //    Math.Pow(((x[6] + _thetaMin * x[7] + Math.Pow(_thetaMin, 2) * x[8] + Math.Pow(_thetaMin, 3) * x[9] + Math.Pow(_thetaMin, 4) * x[10] + Math.Pow(_thetaMin, 5) * x[11]) - (x[6] + _thetaMax * x[7] + Math.Pow(_thetaMax, 2) * x[8] + Math.Pow(_thetaMax, 3) * x[9] + Math.Pow(_thetaMax, 4) * x[10] + Math.Pow(_thetaMax, 5) * x[11])), 2) +
            //    Math.Pow(((x[12] + _thetaMin * x[13] + Math.Pow(_thetaMin, 2) * x[14] + Math.Pow(_thetaMin, 3) * x[15] + Math.Pow(_thetaMin, 4) * x[16] + Math.Pow(_thetaMin, 5) * x[17]) - (x[0] + _thetaMax * x[1] + Math.Pow(_thetaMax, 2) * x[2] + Math.Pow(_thetaMax, 3) * x[3] + Math.Pow(_thetaMax, 4) * x[4] + Math.Pow(_thetaMax, 5) * x[5])), 2);
            func = x[0];
        }



        public void run()
        {
            using (var solver = new NLoptSolver(NLoptAlgorithm.LD_SLSQP, 18, 0.0001, 100))
            {
                solver.SetLowerBounds(new[] { -0.2, -100.0, -100.0, -100.0, -100.0, -100.0, -0.2, -100.0, -100.0, -100.0, -100.0, -100.0, -0.2, -100.0, -100.0, -100.0, -100.0, -100.0 });
                solver.SetUpperBounds(new[] { 0.2, 100.0, 100.0, 100.0, 100.0, 100.0, 0.2, 100.0, 100.0, 100.0, 100.0, 100.0, 0.2, 100.0, 100.0, 100.0, 100.0, 100.0 });

                solver.SetMinObjective(objfun);
                solver.AddLessOrEqualZeroConstraint(alphaMaxConstraint, 0.001);
                solver.AddLessOrEqualZeroConstraint(dthetaConstraint, 0.00001);
                solver.AddEqualZeroConstraint(impactConstraint1, 0.001);
                solver.AddEqualZeroConstraint(impactConstraint2, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint1, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint2, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint3, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint4, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint5, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint6, 0.001);

                double? finalScore;

                var initialValue = _currentParameterValues;
                var result = solver.Optimize(initialValue, out finalScore);
                int len = _currentParameterValues.Length / 3;

                for (int i = 0; i < len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi1Parameters["P" + i.ToString()].RealValue;

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi2Parameters["P" + i.ToString()].RealValue;

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi3Parameters["P" + i.ToString()].RealValue ;

                }

            }

        }
        public void runAlpha()
        {
            using (var solver = new NLoptSolver(NLoptAlgorithm.LD_SLSQP, 18, 0.0001, 100))
            {
                solver.SetLowerBounds(new[] { -0.2, -100.0, -100.0, -100.0, -100.0, -100.0, -0.2, -100.0, -100.0, -100.0, -100.0, -100.0, -0.2, -100.0, -100.0, -100.0, -100.0, -100.0 });
                solver.SetUpperBounds(new[] { 0.2, 100.0, 100.0, 100.0, 100.0, 100.0, 0.2, 100.0, 100.0, 100.0, 100.0, 100.0, 0.2, 100.0, 100.0, 100.0, 100.0, 100.0 });

                solver.SetMinObjective(objfun);
                solver.AddEqualZeroConstraint(geometryConstraint1, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint2, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint3, 0.001);
                //solver.AddEqualZeroConstraint(geometryConstraint4, 0.001);
                //solver.AddEqualZeroConstraint(geometryConstraint5, 0.001);
                //solver.AddEqualZeroConstraint(geometryConstraint6, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint7, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint8, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint9, 0.001);

                double? finalScore;

                var initialValue = _currentParameterValues;
                var result = solver.Optimize(initialValue, out finalScore);
                int len = _currentParameterValues.Length / 3;

                for (int i = 0; i < len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi1Parameters["P" + i.ToString()].RealValue;

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi2Parameters["P" + i.ToString()].RealValue;

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi3Parameters["P" + i.ToString()].RealValue;

                }
            }

        }
        public void runDtheta()
        {
            using (var solver = new NLoptSolver(NLoptAlgorithm.LD_SLSQP, 18, 0.0001, 100))
            {
                solver.SetLowerBounds(new[] { -0.2, -100.0, -100.0, -100.0, -100.0, -100.0, -0.2, -100.0, -100.0, -100.0, -100.0, -100.0, -0.2, -100.0, -100.0, -100.0, -100.0, -100.0 });
                solver.SetUpperBounds(new[] { 0.2, 100.0, 100.0, 100.0, 100.0, 100.0, 0.2, 100.0, 100.0, 100.0, 100.0, 100.0, 0.2, 100.0, 100.0, 100.0, 100.0, 100.0 });

                solver.SetMinObjective(objfun);
                solver.AddLessOrEqualZeroConstraint(alphaMaxConstraint, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint1, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint2, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint3, 0.001);
                //solver.AddEqualZeroConstraint(geometryConstraint4, 0.001);
                //solver.AddEqualZeroConstraint(geometryConstraint5, 0.001);
                //solver.AddEqualZeroConstraint(geometryConstraint6, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint7, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint8, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint9, 0.001);

                double? finalScore;

                var initialValue = _currentParameterValues;
                var result = solver.Optimize(initialValue, out finalScore);
                int len = _currentParameterValues.Length / 3;

                for (int i = 0; i < len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi1Parameters["P" + i.ToString()].RealValue;

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi2Parameters["P" + i.ToString()].RealValue;

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi3Parameters["P" + i.ToString()].RealValue;

                }
            }

        }
        public void runAlphaAndDtheta()
        {
            using (var solver = new NLoptSolver(NLoptAlgorithm.LD_SLSQP, 18, 0.0001, 100))
            {
                solver.SetLowerBounds(new[] { -0.2, -100.0, -100.0, -100.0, -100.0, -100.0, -0.2, -100.0, -100.0, -100.0, -100.0, -100.0, -0.2, -100.0, -100.0, -100.0, -100.0, -100.0 });
                solver.SetUpperBounds(new[] { 0.2, 100.0, 100.0, 100.0, 100.0, 100.0, 0.2, 100.0, 100.0, 100.0, 100.0, 100.0, 0.2, 100.0, 100.0, 100.0, 100.0, 100.0 });

                solver.SetMinObjective(objfunAlphaAndDtheta);
                solver.AddEqualZeroConstraint(geometryConstraint1, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint2, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint3, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint4, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint5, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint6, 0.001);

                double? finalScore;

                var initialValue = _currentParameterValues;
                var result = solver.Optimize(initialValue, out finalScore);
                int len = _currentParameterValues.Length / 3;

                for (int i = 0; i < len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi1Parameters["P" + i.ToString()].RealValue;

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi2Parameters["P" + i.ToString()].RealValue;

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi3Parameters["P" + i.ToString()].RealValue;

                }
            }

        }
        public void runImpact()
        {
            using (var solver = new NLoptSolver(NLoptAlgorithm.LD_SLSQP, 18, 0.0001, 100))
            {
                solver.SetLowerBounds(new[] { -0.2, -100.0, -100.0, -100.0, -100.0, -100.0, -0.2, -100.0, -100.0, -100.0, -100.0, -100.0, -0.2, -100.0, -100.0, -100.0, -100.0, -100.0 });
                solver.SetUpperBounds(new[] { 0.2, 100.0, 100.0, 100.0, 100.0, 100.0, 0.2, 100.0, 100.0, 100.0, 100.0, 100.0, 0.2, 100.0, 100.0, 100.0, 100.0, 100.0 });

                solver.SetMinObjective(objfunImpact);
                solver.AddLessOrEqualZeroConstraint(alphaMaxConstraint, 0.001);
                solver.AddLessOrEqualZeroConstraint(dthetaConstraint, 0.00001);
                solver.AddEqualZeroConstraint(geometryConstraint1, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint2, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint3, 0.001);
                //solver.AddEqualZeroConstraint(geometryConstraint4, 0.001);
                //solver.AddEqualZeroConstraint(geometryConstraint5, 0.001);
                //solver.AddEqualZeroConstraint(geometryConstraint6, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint7, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint8, 0.001);
                solver.AddEqualZeroConstraint(geometryConstraint9, 0.001);

                double? finalScore;

                var initialValue = _currentParameterValues;
                var result = solver.Optimize(initialValue, out finalScore);
                int len = _currentParameterValues.Length / 3;

                for (int i = 0; i < len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi1Parameters["P" + i.ToString()].RealValue;

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi2Parameters["P" + i.ToString()].RealValue;

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _currentParameterValues[i] = _gait.vhc.phi3Parameters["P" + i.ToString()].RealValue;

                }
            }

        }

        public double objfun(double[] p, double[] grad)
        {
            int len = p.Length / 3;

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
            double[] values = evalTorques(_thetaMin, _thetaMax);

            if (grad != null)
            {
                double[] gradValues;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValues = evalTorques(_thetaMin, _thetaMax);
                    grad[i] = (Math.Pow(gradValues[0], 2) + Math.Pow(gradValues[1], 2) - (Math.Pow(values[0], 2) + Math.Pow(values[1], 2))) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValues = evalTorques(_thetaMin, _thetaMax);
                    grad[i] = (Math.Pow(gradValues[0], 2) + Math.Pow(gradValues[1], 2) - (Math.Pow(values[0], 2) + Math.Pow(values[1], 2))) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValues = evalTorques(_thetaMin, _thetaMax);
                    grad[i] = (Math.Pow(gradValues[0], 2) + Math.Pow(gradValues[1], 2) - (Math.Pow(values[0], 2) + Math.Pow(values[1], 2))) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return Math.Pow(values[0], 2) + Math.Pow(values[1], 2);

        }
        public double objfunAlpha(double[] p, double[] grad)
        {
            int len = p.Length / 3;
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
            double alpha = evaluateAlphaConstraint(_thetaMin, _thetaMax);
            if (grad != null)
            {
                double alphaVal;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    alphaVal = evaluateAlphaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (alphaVal - (alpha)) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    alphaVal = evaluateAlphaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (alphaVal - (alpha)) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    alphaVal = evaluateAlphaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (alphaVal - (alpha)) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return alpha;

        }
        public double objfunDtheta(double[] p, double[] grad)
        {
            int len = p.Length / 3;
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
            double dthetaMin = evaluateDthetaConstraint(_thetaMin, _thetaMax);
            if (grad != null)
            {
                double dthetaMinVal;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    dthetaMinVal = evaluateDthetaConstraint(_thetaMin, _thetaMax);
                    grad[i] = ((-1) * dthetaMinVal - ((-1) * dthetaMin)) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    dthetaMinVal = evaluateDthetaConstraint(_thetaMin, _thetaMax);
                    grad[i] = ((-1) * dthetaMinVal - ((-1) * dthetaMin)) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    dthetaMinVal = evaluateDthetaConstraint(_thetaMin, _thetaMax);
                    grad[i] = ((-1) * dthetaMinVal - ((-1) * dthetaMin)) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return (-1) * dthetaMin;

        }
        public double objfunAlphaAndDtheta(double[] p, double[] grad)
        {
            int len = p.Length / 3;
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
            double alpha = evaluateAlphaConstraint(_thetaMin, _thetaMax);
            double dthetaMin = evaluateDthetaConstraint(_thetaMin, _thetaMax);
            if (grad != null)
            {
                double alphaVal;
                double dthetaMinVal;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    alphaVal = evaluateAlphaConstraint(_thetaMin, _thetaMax);
                    dthetaMinVal = evaluateDthetaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (alphaVal + (-1) * dthetaMinVal - (alpha + (-1) * dthetaMin)) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    alphaVal = evaluateAlphaConstraint(_thetaMin, _thetaMax);
                    dthetaMinVal = evaluateDthetaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (alphaVal + (-1) * dthetaMinVal - (alpha + (-1) * dthetaMin)) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    alphaVal = evaluateAlphaConstraint(_thetaMin, _thetaMax);
                    dthetaMinVal = evaluateDthetaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (alphaVal + (-1) * dthetaMinVal - (alpha + (-1) * dthetaMin)) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return alpha + (-1) * dthetaMin;

        }
        public double objfunImpact(double[] p, double[] grad)
        {
            int len = p.Length / 3;

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
            double impactValue1 = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactSecondLine(_thetaMin, _thetaMax);
            double impactValue2 = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactThirdLine(_thetaMin, _thetaMax);
            if (grad != null)
            {
                double gradValue1;
                double gradValue2;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue1 = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactSecondLine(_thetaMin, _thetaMax);
                    gradValue2 = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactThirdLine(_thetaMin, _thetaMax);
                    grad[i] = (Math.Pow(gradValue1, 2) + Math.Pow(gradValue2, 2) - (Math.Pow(impactValue1, 2) + Math.Pow(impactValue2, 2))) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue1 = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactSecondLine(_thetaMin, _thetaMax);
                    gradValue2 = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactThirdLine(_thetaMin, _thetaMax);
                    grad[i] = (Math.Pow(gradValue1, 2) + Math.Pow(gradValue2, 2) - (Math.Pow(impactValue1, 2) + Math.Pow(impactValue2, 2))) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue1 = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactSecondLine(_thetaMin, _thetaMax);
                    gradValue2 = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactThirdLine(_thetaMin, _thetaMax);
                    grad[i] = (Math.Pow(gradValue1, 2) + Math.Pow(gradValue2, 2) - (Math.Pow(impactValue1, 2) + Math.Pow(impactValue2, 2))) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return Math.Pow(impactValue1, 2) + Math.Pow(impactValue2, 2);

        }


        public double geometryConstraint1(double[] p, double[] grad)
        {

            if (grad != null)
            {
                grad[0] = 1;
                grad[1] = _thetaMin;
                grad[2] = Math.Pow(_thetaMin, 2);
                grad[3] = Math.Pow(_thetaMin, 3);
                grad[4] = Math.Pow(_thetaMin, 4);
                grad[5] = Math.Pow(_thetaMin, 5);
                grad[6] = 0;
                grad[7] = 0;
                grad[8] = 0;
                grad[9] = 0;
                grad[10] = 0;
                grad[11] = 0;
                grad[12] = -1;
                grad[13] = -_thetaMax;
                grad[14] = -Math.Pow(_thetaMax, 2);
                grad[15] = -Math.Pow(_thetaMax, 3);
                grad[16] = -Math.Pow(_thetaMax, 4);
                grad[17] = -Math.Pow(_thetaMax, 5);
            }
            return p[0] + p[1] * _thetaMin + p[2] * Math.Pow(_thetaMin, 2) + p[3] * Math.Pow(_thetaMin, 3) + p[4] * Math.Pow(_thetaMin, 4) + p[5] * Math.Pow(_thetaMin, 5) - 
                (p[12] + p[13] * _thetaMax + p[14] * Math.Pow(_thetaMax, 2) + p[15] * Math.Pow(_thetaMax, 3) + p[16] * Math.Pow(_thetaMax, 4) + p[17] * Math.Pow(_thetaMax, 5));
        }
        public double geometryConstraint2(double[] p, double[] grad)
        {
            if (grad != null)
            {
                grad[0] = 0;
                grad[1] = 0;
                grad[2] = 0;
                grad[3] = 0;
                grad[4] = 0;
                grad[5] = 0;
                grad[6] = 0;
                grad[7] = _thetaMin - _thetaMax;
                grad[8] = Math.Pow(_thetaMin, 2) - Math.Pow(_thetaMax, 2);
                grad[9] = Math.Pow(_thetaMin, 3) - Math.Pow(_thetaMax, 3);
                grad[10] = Math.Pow(_thetaMin, 4) - Math.Pow(_thetaMax, 4);
                grad[11] = Math.Pow(_thetaMin, 5) - Math.Pow(_thetaMax, 5);
                grad[12] = 0;
                grad[13] = 0;
                grad[14] = 0;
                grad[15] = 0;
                grad[16] = 0;
                grad[17] = 0;
            }
            return p[6] + p[7] * _thetaMin + p[8] * Math.Pow(_thetaMin, 2) + p[9] * Math.Pow(_thetaMin, 3) + p[10] * Math.Pow(_thetaMin, 4) + p[11] * Math.Pow(_thetaMin, 5) -
                (p[6] + p[7] * _thetaMax + p[8] * Math.Pow(_thetaMax, 2) + p[9] * Math.Pow(_thetaMax, 3) + p[10] * Math.Pow(_thetaMax, 4) + p[11] * Math.Pow(_thetaMax, 5));
        }
        public double geometryConstraint3(double[] p, double[] grad)
        {
            if (grad != null)
            {
                grad[0] = -1;
                grad[1] = -_thetaMax;
                grad[2] = -Math.Pow(_thetaMax, 2);
                grad[3] = -Math.Pow(_thetaMax, 3);
                grad[4] = -Math.Pow(_thetaMax, 4);
                grad[5] = -Math.Pow(_thetaMax, 5);
                grad[6] = 0;
                grad[7] = 0;
                grad[8] = 0;
                grad[9] = 0;
                grad[10] = 0;
                grad[11] = 0;
                grad[12] = 1;
                grad[13] = _thetaMin;
                grad[14] = Math.Pow(_thetaMin, 2);
                grad[15] = Math.Pow(_thetaMin, 3);
                grad[16] = Math.Pow(_thetaMin, 4);
                grad[17] = Math.Pow(_thetaMin, 5);
            }
            return p[12] + p[13] * _thetaMin + p[14] * Math.Pow(_thetaMin, 2) + p[15] * Math.Pow(_thetaMin, 3) + p[16] * Math.Pow(_thetaMin, 4) + p[17] * Math.Pow(_thetaMin, 5) -
                (p[0] + p[1] * _thetaMax + p[2] * Math.Pow(_thetaMax, 2) + p[3] * Math.Pow(_thetaMax, 3) + p[4] * Math.Pow(_thetaMax, 4) + p[5] * Math.Pow(_thetaMax, 5));
        }

        public double geometryConstraint4(double[] p, double[] grad)
        {

            if (grad != null)
            {
                grad[0] = 1 + 1;
                grad[1] = _thetaMin + _thetaMax;
                grad[2] = Math.Pow(_thetaMin, 2) + Math.Pow(_thetaMax, 2);
                grad[3] = Math.Pow(_thetaMin, 3) + Math.Pow(_thetaMax, 3);
                grad[4] = Math.Pow(_thetaMin, 4) + Math.Pow(_thetaMax, 4);
                grad[5] = Math.Pow(_thetaMin, 5) + Math.Pow(_thetaMax, 5);
                grad[6] = 0;
                grad[7] = 0;
                grad[8] = 0;
                grad[9] = 0;
                grad[10] = 0;
                grad[11] = 0;
                grad[12] = 0;
                grad[13] = 0;
                grad[14] = 0;
                grad[15] = 0;
                grad[16] = 0;
                grad[17] = 0;
            }
            return p[0] + p[1] * _thetaMin + p[2] * Math.Pow(_thetaMin, 2) + p[3] * Math.Pow(_thetaMin, 3) + p[4] * Math.Pow(_thetaMin, 4) + p[5] * Math.Pow(_thetaMin, 5) +
                p[0] + p[1] * _thetaMax + p[2] * Math.Pow(_thetaMax, 2) + p[3] * Math.Pow(_thetaMax, 3) + p[4] * Math.Pow(_thetaMax, 4) + p[5] * Math.Pow(_thetaMax, 5);
        }
        public double geometryConstraint5(double[] p, double[] grad)
        {
            if (grad != null)
            {
                grad[0] = 0;
                grad[1] = 0;
                grad[2] = 0;
                grad[3] = 0;
                grad[4] = 0;
                grad[5] = 0;
                grad[6] = 0;
                grad[7] = _thetaMin - _thetaMax;
                grad[8] = Math.Pow(_thetaMin, 2) - Math.Pow(_thetaMax, 2);
                grad[9] = Math.Pow(_thetaMin, 3) - Math.Pow(_thetaMax, 3);
                grad[10] = Math.Pow(_thetaMin, 4) - Math.Pow(_thetaMax, 4);
                grad[11] = Math.Pow(_thetaMin, 5) - Math.Pow(_thetaMax, 5);
                grad[12] = 0;
                grad[13] = 0;
                grad[14] = 0;
                grad[15] = 0;
                grad[16] = 0;
                grad[17] = 0;
            }
            return p[6] + p[7] * _thetaMin + p[8] * Math.Pow(_thetaMin, 2) + p[9] * Math.Pow(_thetaMin, 3) + p[10] * Math.Pow(_thetaMin, 4) + p[11] * Math.Pow(_thetaMin, 5) -
                p[6] + p[7] * _thetaMax + p[8] * Math.Pow(_thetaMax, 2) + p[9] * Math.Pow(_thetaMax, 3) + p[10] * Math.Pow(_thetaMax, 4) + p[11] * Math.Pow(_thetaMax, 5);
        }
        public double geometryConstraint6(double[] p, double[] grad)
        {
            if (grad != null)
            {
                grad[0] = 0;
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
                grad[12] = 1 + 1;
                grad[13] = _thetaMin + _thetaMax;
                grad[14] = Math.Pow(_thetaMin, 2) + Math.Pow(_thetaMax, 2);
                grad[15] = Math.Pow(_thetaMin, 3) + Math.Pow(_thetaMax, 3);
                grad[16] = Math.Pow(_thetaMin, 4) + Math.Pow(_thetaMax, 4);
                grad[17] = Math.Pow(_thetaMin, 5) + Math.Pow(_thetaMax, 5);
            }
            return p[12] + p[13] * _thetaMin + p[14] * Math.Pow(_thetaMin, 2) + p[15] * Math.Pow(_thetaMin, 3) + p[16] * Math.Pow(_thetaMin, 4) + p[17] * Math.Pow(_thetaMin, 5) +
                p[12] + p[13] * _thetaMax + p[14] * Math.Pow(_thetaMax, 2) + p[15] * Math.Pow(_thetaMax, 3) + p[16] * Math.Pow(_thetaMax, 4) + p[17] * Math.Pow(_thetaMax, 5);
        }

        public double geometryConstraint7(double[] p, double[] grad)
        {

            if (grad != null)
            {
                grad[0] = 1;
                grad[1] = _thetaMin;
                grad[2] = Math.Pow(_thetaMin, 2);
                grad[3] = Math.Pow(_thetaMin, 3);
                grad[4] = Math.Pow(_thetaMin, 4);
                grad[5] = Math.Pow(_thetaMin, 5);
                grad[6] = 0;
                grad[7] = 0;
                grad[8] = 0;
                grad[9] = 0;
                grad[10] = 0;
                grad[11] = 0;
                grad[12] = 0;
                grad[13] = 0;
                grad[14] = 0;
                grad[15] = 0;
                grad[16] = 0;
                grad[17] = 0;
            }
            return p[0] + p[1] * _thetaMin + p[2] * Math.Pow(_thetaMin, 2) + p[3] * Math.Pow(_thetaMin, 3) + p[4] * Math.Pow(_thetaMin, 4) + p[5] * Math.Pow(_thetaMin, 5) - _q1;
        }
        public double geometryConstraint8(double[] p, double[] grad)
        {
            if (grad != null)
            {
                grad[0] = 0;
                grad[1] = 0;
                grad[2] = 0;
                grad[3] = 0;
                grad[4] = 0;
                grad[5] = 0;
                grad[6] = 1;
                grad[7] = _thetaMin;
                grad[8] = Math.Pow(_thetaMin, 2);
                grad[9] = Math.Pow(_thetaMin, 3);
                grad[10] = Math.Pow(_thetaMin, 4);
                grad[11] = Math.Pow(_thetaMin, 5);
                grad[12] = 0;
                grad[13] = 0;
                grad[14] = 0;
                grad[15] = 0;
                grad[16] = 0;
                grad[17] = 0;
            }
            return p[6] + p[7] * _thetaMin + p[8] * Math.Pow(_thetaMin, 2) + p[9] * Math.Pow(_thetaMin, 3) + p[10] * Math.Pow(_thetaMin, 4) + p[11] * Math.Pow(_thetaMin, 5) - _q2;
        }
        public double geometryConstraint9(double[] p, double[] grad)
        {
            if (grad != null)
            {
                grad[0] = 0;
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
                grad[12] = 1;
                grad[13] = _thetaMin;
                grad[14] = Math.Pow(_thetaMin, 2);
                grad[15] = Math.Pow(_thetaMin, 3);
                grad[16] = Math.Pow(_thetaMin, 4);
                grad[17] = Math.Pow(_thetaMin, 5);
            }
            return p[12] + p[13] * _thetaMin + p[14] * Math.Pow(_thetaMin, 2) + p[15] * Math.Pow(_thetaMin, 3) + p[16] * Math.Pow(_thetaMin, 4) + p[17] * Math.Pow(_thetaMin, 5) - _q3;
        }

        public double alphaMaxConstraint(double[] p, double[] grad)
        {
            int len = p.Length / 3;

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
            double alpha = evaluateAlphaConstraint(_thetaMin, _thetaMax);
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateAlphaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (gradValue - alpha) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateAlphaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (gradValue - alpha) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateAlphaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (gradValue - alpha) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return alpha;
        }
        public double dthetaConstraint(double[] p, double[] grad)
        {
            int len = p.Length / 3;

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
            double dthetaMin = evaluateDthetaConstraint(_thetaMin, _thetaMax);
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateDthetaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (gradValue - dthetaMin) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateDthetaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (gradValue - dthetaMin) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateDthetaConstraint(_thetaMin, _thetaMax);
                    grad[i] = (gradValue - dthetaMin) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return (-1) * dthetaMin;
        }
        public double impactConstraint1(double[] p, double[] grad)
        {
            int len = p.Length / 3;

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
            double impactValue = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactSecondLine(_thetaMin, _thetaMax);
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactSecondLine(_thetaMin, _thetaMax);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactSecondLine(_thetaMin, _thetaMax);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactSecondLine(_thetaMin, _thetaMax);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return impactValue;
        }
        public double impactConstraint2(double[] p, double[] grad)
        {
            int len = p.Length / 3;

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
            double impactValue = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactThirdLine(_thetaMin, _thetaMax);
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactThirdLine(_thetaMin, _thetaMax);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactThirdLine(_thetaMin, _thetaMax);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(_thetaMin, _thetaMax) - _gait.impactThirdLine(_thetaMin, _thetaMax);
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

        public void thetaRange(double q1, double q2, double q3)
        {
            double a = -(_param.m1 * _param.l1 + _param.m2 * _param.L1 + _param.m3 * _param.L1) / (_param.m1 + _param.m2 + _param.m3);
            double b = -(_param.m2 * _param.l2) / (_param.m1 + _param.m2 + _param.m3);
            double c = (_param.m3 * _param.l3) / (_param.m1 + _param.m2 + _param.m3);
            _thetaMin = a * q1 + b * q2 + c * q3;
            _thetaMax = a * (-q1) + b * (q2) + c * (-q3);
        }
        public double evaluateAlphaConstraint(double thetaMin, double thetaMax)
        {

            double steps = 50;
            double dx = (thetaMax - thetaMin) / steps;
            double theta = thetaMin;
            double val = double.MinValue;
            double alphaVal = 0;
            for (int i = 0; i < steps; i++)
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
