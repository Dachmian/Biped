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

        double _q1End;
        double _q2End;
        double _q3End;

        private double[] _currentParameterValues;
        private static double _step = 1e-9;

        public delegate double func(Expression exp);

        public SLSQPCoM(BRgait gait, int numOfParams)
        {
            _gait = gait;
            BRVHC vhc = gait.vhc;
            _currentParameterValues = new double[3 * numOfParams];
        }







        public void runNumericalSLSQP()
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
            double[] values = evalTorques();

            if (grad != null)
            {
                double[] gradValues;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValues = evalTorques();
                    grad[i] = (Math.Pow(gradValues[0], 2) + Math.Pow(gradValues[1], 2) - (Math.Pow(values[0], 2) + Math.Pow(values[1], 2))) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValues = evalTorques();
                    grad[i] = (Math.Pow(gradValues[0], 2) + Math.Pow(gradValues[1], 2) - (Math.Pow(values[0], 2) + Math.Pow(values[1], 2))) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValues = evalTorques();
                    grad[i] = (Math.Pow(gradValues[0], 2) + Math.Pow(gradValues[1], 2) - (Math.Pow(values[0], 2) + Math.Pow(values[1], 2))) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return Math.Pow(values[0], 2) + Math.Pow(values[1], 2);

        }
        public double geometryConstraint1(double[] p, double[] grad)
        {
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
                grad[12] = -1;
                grad[13] = -1;
                grad[14] = -1;
                grad[15] = -1;
                grad[16] = -1;
                grad[17] = -1;
            }
            return p[0] - (p[12] + p[13] + p[14] + p[15] + p[16] + p[17]);
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
                grad[7] = -1;
                grad[8] = -1;
                grad[9] = -1;
                grad[10] = -1;
                grad[11] = -1;
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
            if (grad != null)
            {
                grad[0] = -1;
                grad[1] = -1;
                grad[2] = -1;
                grad[3] = -1;
                grad[4] = -1;
                grad[5] = -1;
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
            double alpha = evaluateAlphaConstraint();
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateAlphaConstraint();
                    grad[i] = (gradValue - alpha) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateAlphaConstraint();
                    grad[i] = (gradValue - alpha) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateAlphaConstraint();
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
            double dthetaMin = evaluateDthetaConstraint();
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateDthetaConstraint();
                    grad[i] = (gradValue - dthetaMin) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateDthetaConstraint();
                    grad[i] = (gradValue - dthetaMin) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = evaluateDthetaConstraint();
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
            double impactValue = _gait.impactFirstLine(0, 1) - _gait.impactSecondLine(0, 1);
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(0, 1) - _gait.impactSecondLine(0, 1);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(0, 1) - _gait.impactSecondLine(0, 1);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(0, 1) - _gait.impactSecondLine(0, 1);
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
            double impactValue = _gait.impactFirstLine(0, 1) - _gait.impactThirdLine(0, 1);
            if (grad != null)
            {
                double gradValue;
                for (int i = 0; i < len; i++)
                {
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(0, 1) - _gait.impactThirdLine(0, 1);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi1Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = len; i < 2 * len; i++)
                {
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(0, 1) - _gait.impactThirdLine(0, 1);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi2Parameters["P" + i.ToString()] = p[i];

                }
                for (int i = 2 * len; i < 3 * len; i++)
                {
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i] + _step;
                    gradValue = _gait.impactFirstLine(0, 1) - _gait.impactThirdLine(0, 1);
                    grad[i] = (gradValue - impactValue) / _step;
                    _gait.vhc.phi3Parameters["P" + i.ToString()] = p[i];

                }
            }
            return impactValue;

        }
        public double[] evalTorques()
        {
            double[][] firstIntegral = TrapezoidalSum2.calculateFirstIntegral(_gait.vhc.evalTwoTimesBetaDividedByAlpha, 0, 1, 200);
            double[][] secondIntegral = TrapezoidalSum2.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1], 0, 1, 200);
            int len = firstIntegral[1].Length;
            double dthetaTSquared = (-secondIntegral[1][len - 1] * Math.Exp(firstIntegral[1][len - 1])) / (1 - Math.Exp(firstIntegral[1][len - 1]) * Math.Pow(_gait.impactSecondLine(0, 1), 2));
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(0, 1), 2);
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
                THETA[0, 0] = 0;
                THETA[1, 0] = dtheta0Squared;

                double torque1 = 0;
                double torque2 = 0;
                double theta = 0;
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
                double ddtheta0 = (-_gait.vhc.evalBeta(0) * dtheta0Squared - _gait.vhc.evalGamma(0)) / _gait.vhc.evalAlpha(0);
                torque1Diff = (_gait.vhc.evalAlpha1(0) * ddtheta0 + _gait.vhc.evalBeta1(0) * dtheta0Squared + _gait.vhc.evalGamma1(0)) -
                    (_gait.vhc.evalAlpha3(1) * ddtheta + _gait.vhc.evalBeta3(1) * dthetaTSquared + _gait.vhc.evalGamma3(1));

                torque2Diff = _gait.vhc.evalAlpha3(0) * ddtheta0 + _gait.vhc.evalBeta3(0) * dtheta0Squared + _gait.vhc.evalGamma3(0) -
                    (_gait.vhc.evalAlpha1(1) * ddtheta + _gait.vhc.evalBeta1(1) * dthetaTSquared + _gait.vhc.evalGamma1(1));
                values = new double[] { torque1Max, torque2Max, torque1Diff, torque2Diff };

            }
            return values;
        }

        public double evaluateAlphaConstraint()
        {


            double dx = 0.02;
            double theta = 0;
            double val = double.MinValue;
            double alphaVal = 0;
            for (int i = 0; i < 1 / dx; i++)
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
        public double evaluateDthetaConstraint()
        {
            double[][] firstIntegral = TrapezoidalSum2.calculateFirstIntegral(_gait.vhc.evalTwoTimesBetaDividedByAlpha, 0, 1, 200);
            double[][] secondIntegral = TrapezoidalSum2.calculateSecondIntegral(_gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral[1], 0, 1, 200);
            int len = firstIntegral[1].Length;
            double dthetaTSquared = (-secondIntegral[1][len - 1] * Math.Exp(firstIntegral[1][len - 1])) / (1 - Math.Exp(firstIntegral[1][len - 1]) * Math.Pow(_gait.impactSecondLine(0, 1), 2));
            double dtheta0Squared = dthetaTSquared * Math.Pow(_gait.impactSecondLine(0, 1), 2);
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
                THETA[0, 0] = 0;
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
