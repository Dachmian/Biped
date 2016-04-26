using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Symbolics;
using MathNet.Numerics.LinearAlgebra;
using System.IO;
using System.Globalization;

namespace BipedRobot{
    static class Program{
        static void Main(string[] args)
        {
            Biped biped = new Biped (@"../../basic_randomized_set.xml");
            //IntegrationFullDynamics.run(ref biped);
            //plotting.plotStates(biped);
            //gaitSearch.run(ref biped);
            //findGaitsManually(biped);
            findGaitsHere(biped);
            
        }

        static void test()
        {
            Dictionary<string, FloatingPoint>  _parameters = new Dictionary<string, FloatingPoint>();
            _parameters.Add("g", 1);
            _parameters.Add("m1", 1);
            _parameters.Add("m2", 1);
            _parameters.Add("m3", 1);

            _parameters.Add("l1", 1);
            _parameters.Add("l2", 1);
            _parameters.Add("l3", 10);

            _parameters.Add("L1", 10);
            _parameters.Add("L2", 1);
            _parameters.Add("L3", 1);

            _parameters.Add("J1", 10);
            _parameters.Add("J2", 1);
            _parameters.Add("J3", 10);
            _parameters.Add("theta", 0);
            _parameters.Add("dtheta", 0);
            _parameters.Add("ddtheta", 0);
            _parameters.Add("P0", 1);
            _parameters.Add("P1", 1);
            _parameters.Add("P2", 0);
            _parameters.Add("P3", 0);
            _parameters.Add("P4", 0);
            _parameters.Add("P5", 0);
            _parameters.Add("P6", 0);
            _parameters.Add("P7", 0);
            _parameters.Add("P8", 0);
            _parameters.Add("P9", 0);
            _parameters.Add("P10", 0);
            _parameters.Add("P11", 0);
            _parameters.Add("P12", 0);
            _parameters.Add("P13", 0);
            _parameters.Add("P14", 1);

            Expression test = Infix.ParseOrUndefined("2*(-((J1 - L1*l1*m1 + l1^2*m1)*(0^4*L1*L3*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) - 0^4*L3*l1*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)))/((4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1) - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4))^2 + (1^4*L1*l2*m1*sin(1^4*P0 + 1^3*P1 + 1^2*P2 + 1*P3 + P4 - (1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)))/((4*0^3*P5 + 3*0^2*P6 + 2*0*P7 + P8)*(J2 + l2^2*m1) + L3*l2*m1*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9)))*((J1 - L1*l1*m1 + l1^2*m1)/((4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1) - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)) - ((4*1^3*P5 + 3*1^2*P6 + 2*1*P7 + P8)*(J2 + l2^2*m1) + L1*l2*m1*cos(1^4*P0 + 1^3*P1 + 1^2*P2 + 1*P3 + P4 - (1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)))/((4*0^3*P5 + 3*0^2*P6 + 2*0*P7 + P8)*(J2 + l2^2*m1) + L3*l2*m1*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9)))^(-1 + 2) - 2*((J1 - L1*l1*m1 + l1^2*m1)/((4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1) - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)) - (J1 - L1*l1*m1 + l1^2*m1 + (4*1^3*P10 + 3*1^2*P11 + 2*1*P12 + P13)*(J3 - L3*l3*m3 + l3^2*m3) + L3*l1*m1*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + L1*L3*m2*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + L1*L3*m3*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) - L1*l3*m3*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + (4*1^3*P5 + 3*1^2*P6 + 2*1*P7 + P8)*(J2 + l2^2*m2 + L3*l2*m2*cos(-(1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)) + L1*l2*m2*cos(1^4*P0 + 1^3*P1 + 1^2*P2 + 1*P3 + P4 - (1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)))/(J3 + L3^2*m1 + L3^2*m2 + L3^2*m3 - 2*L3*l3*m3 + l3^2*m3 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + (4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)) + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9) + (4*0^3*P5 + 3*0^2*P6 + 2*0*P7 + P8)*(J2 + l2^2*m2 + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9))))^(-1 + 2)*(-((J1 - L1*l1*m1 + l1^2*m1)*(0^4*L1*L3*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) - 0^4*L3*l1*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)))/((4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1) - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4))^2 - ((-1^4*L3*l1*m1*sin(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) - 1^4*L1*L3*m2*sin(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) - 1^4*L1*L3*m3*sin(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + 1^4*L1*l3*m3*sin(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) - 1^4*L1*l2*m2*sin(1^4*P0 + 1^3*P1 + 1^2*P2 + 1*P3 + P4 - (1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)))/(J3 + L3^2*m1 + L3^2*m2 + L3^2*m3 - 2*L3*l3*m3 + l3^2*m3 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + (4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)) + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9) + (4*0^3*P5 + 3*0^2*P6 + 2*0*P7 + P8)*(J2 + l2^2*m2 + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9))) - ((0^4*L1*L3*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) - 0^4*L3*l1*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + (4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(0^4*L1*L3*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) - 0^4*L3*l1*m1*sin(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)))*(J1 - L1*l1*m1 + l1^2*m1 + (4*1^3*P10 + 3*1^2*P11 + 2*1*P12 + P13)*(J3 - L3*l3*m3 + l3^2*m3) + L3*l1*m1*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + L1*L3*m2*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + L1*L3*m3*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) - L1*l3*m3*cos(1^4*P0 + 1^3*P1 - (1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^2*P2 + 1*P3 + P4) + (4*1^3*P5 + 3*1^2*P6 + 2*1*P7 + P8)*(J2 + l2^2*m2 + L3*l2*m2*cos(-(1^4*P10 + 1^3*P11 + 1^2*P12 + 1*P13 + P14) + 1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9)) + L1*l2*m2*cos(1^4*P0 + 1^3*P1 + 1^2*P2 + 1*P3 + P4 - (1^4*P5 + 1^3*P6 + 1^2*P7 + 1*P8 + P9))))/(J3 + L3^2*m1 + L3^2*m2 + L3^2*m3 - 2*L3*l3*m3 + l3^2*m3 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + (4*0^3*P10 + 3*0^2*P11 + 2*0*P12 + P13)*(J1 + L1^2*m1 - 2*L1*l1*m1 + l1^2*m1 - L1*L3*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4) + L3*l1*m1*cos(0^4*P0 + 0^3*P1 - (0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^2*P2 + 0*P3 + P4)) + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9) + (4*0^3*P5 + 3*0^2*P6 + 2*0*P7 + P8)*(J2 + l2^2*m2 + L3*l2*m2*cos(-(0^4*P10 + 0^3*P11 + 0^2*P12 + 0*P13 + P14) + 0^4*P5 + 0^3*P6 + 0^2*P7 + 0*P8 + P9)))^2))");
            Console.WriteLine(MathNet.Symbolics.Evaluate.Evaluate(_parameters, test).RealValue);
            
        }

        static void findGaitsManually(Biped biped)
        {
            int numberOfPoints = 4;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);


            string[] parameters = File.ReadAllLines(@"../../../parameters.txt");
            double dtheta0 = Convert.ToDouble(File.ReadAllLines(@"../../../dtheta0.txt")[0], CultureInfo.InvariantCulture);
            double dthetaT = Convert.ToDouble(File.ReadAllLines(@"../../../dthetaT.txt")[0], CultureInfo.InvariantCulture);
            for (int i = 0; i < parameters.Length; i++)
            {
                Console.WriteLine(parameters[i]);
            }
            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);

            for (int i = 0; i < 4; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi1Parameters = tuple.Item1;

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 4; i < 8; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi2Parameters = tuple.Item1;

            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 8; i < 12; i++)
            {
                tuple.Item1["P" + i.ToString()] = Convert.ToDouble(parameters[i], CultureInfo.InvariantCulture);

            }
            gait.vhc.phi3Parameters = tuple.Item1;

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined(brCrv.dphi2ToString());
            gait.vhc.dphi3 = Infix.ParseOrUndefined(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = Infix.ParseOrUndefined(brCrv.ddphi3ToString());




            Console.WriteLine(gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactSecondLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactThirdLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));

            Console.WriteLine(gait.vhc.evalAlpha(0));
            Console.WriteLine(gait.vhc.evalAlpha(0.5));
            Console.WriteLine(gait.vhc.evalAlpha(1));


            double[,] THETA = evalDthetaConstraint(gait, Math.Pow(dtheta0, 2), Math.Pow(dthetaT, 2));

            Phaseportrait plot = new Phaseportrait(THETA);

            BRReducedData data = integrationReducedDynamics.run(gait.vhc, Vector<double>.Build.Dense(new double[] { 0.0, dtheta0 }),
                Vector<double>.Build.Dense(new double[] { 1, dthetaT }));
            biped.reducedData = data;

            graph graph = new graph(biped);
        }

        public static double[,] evalDthetaConstraint(BRgait gait, double dtheta0Squared, double dthetaTSquared)
        {
            double[,] firstIntegral = RiemannSum.calculateFirstIntegral(gait.vhc.evalTwoTimesBetaDividedByAlpha);
            double[,] secondIntegral = RiemannSum.calculateSecondIntegral(gait.vhc.evalTwoTimesGammaDividedByAlpha, firstIntegral);

            double endi1 = firstIntegral[499, 1];
            double endi2 = secondIntegral[499, 1];
            int len = secondIntegral.Length / secondIntegral.Rank;


            
            double[,] THETA = new double[len,2];
            THETA[0, 0] = 0;
            THETA[0, 1] = dtheta0Squared;
            for (int i = 1; i < len; i++)
            {
                THETA[i, 0] = secondIntegral[i, 0];
                THETA[i, 1] = -secondIntegral[i, 1] + Math.Exp(-firstIntegral[i, 1]) * dtheta0Squared;
            }
            return THETA;
        }



        static void findGaitsHere(Biped biped)
        {

            int numberOfPoints = 4;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);

            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);

            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined(brCrv.dphi2ToString());
            gait.vhc.dphi3 = Infix.ParseOrUndefined(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = Infix.ParseOrUndefined(brCrv.ddphi3ToString());
            double m1 = biped.param.m1;
            double m2 = biped.param.m2;
            double m3 = biped.param.m3;

            double l1 = biped.param.l1;
            double l2 = biped.param.l2;
            double l3 = biped.param.l3;

            double L1 = biped.param.L1;
            double L2 = biped.param.L2;
            double L3 = biped.param.L3;

            double J1 = biped.param.J1;
            double J2 = biped.param.J2;
            double J3 = biped.param.J3;

            double startq1 = -Math.PI / 26;
            double startq2 = -Math.PI / 34;
            double startq3 = Math.PI / 26;

            double endq1 = -startq1;
            double endq2 = startq2;
            double endq3 = -startq3;

            double pLHS00 = L1 * L1 * m1 - 2 * L1 * l1 * m1 + m1 * l1 * l1 + J1;
            double pLHS01 = 0;
            double pLHS02 = -L1 * L3 * m1 * Math.Cos(endq1 - endq3) + L3 * l1 * m1 * Math.Cos(endq1 - endq3);

            double pLHS10 = 0;
            double pLHS11 = l2 * l2 * m2 + J2;
            double pLHS12 = L3 * l2 * m2 * Math.Cos(endq2 - endq3);

            double pLHS20 = L1 * L1 * m1 - L1 * L3 * m1 * Math.Cos(endq1 - endq3) - 2 * L1 * l1 * m1 + L3 * l1 * m1 * Math.Cos(endq1 - endq3) + m1 * l1 * l1 + J1;
            double pLHS21 = L3 * l2 * m2 * Math.Cos(endq2 - endq3) + l2 * l2 * m2 + J2;
            double pLHS22 = L3 * l1 * m1 * Math.Cos(endq1 - endq3) + L3 * L3 * m1 + L3 * L3 * m2 + L3 * L3 * m3 - L1 * L3 * m1 * Math.Cos(endq1 - endq3) + L3 * l2 * m2 * Math.Cos(endq2 - endq3) - 0.2e1 * m3 * l3 * L3 + l3 * l3 * m3 + J3;


            double pRHS00 = -L1 * l1 * m1 + m1 * l1 * l1 + J1;
            double pRHS01 = 0;
            double pRHS02 = 0;

            double pRHS10 = L1 * Math.Cos(endq1 - endq2) * l2 * m2;
            double pRHS11 = l2 * l2 * m2 + J2;
            double pRHS12 = 0;

            double pRHS20 = J1 + L3 * l1 * m1 * Math.Cos(endq1 - endq3) - L1 * l1 * m1 + L1 * L3 * m2 * Math.Cos(endq1 - endq3) + m3 * L1 * L3 * Math.Cos(endq1 - endq3) - L1 * Math.Cos(endq1 - endq3) * l3 * m3 + L1 * Math.Cos(endq1 - endq2) * l2 * m2 + m1 * l1 * l1;
            double pRHS21 = L3 * l2 * m2 * Math.Cos(endq2 - endq3) + l2 * l2 * m2 + J2;
            double pRHS22 = -m3 * l3 * L3 + l3 * l3 * m3 + J3;

            Matrix<double> pLHS = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {pLHS00,pLHS01,pLHS02 },
                {pLHS10,pLHS11,pLHS12 },
                {pLHS20,pLHS21,pLHS22 },
            });

            Matrix<double> pRHS = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {pRHS00,pRHS01,pRHS02 },
                {pRHS10,pRHS11,pRHS12 },
                {pRHS20,pRHS21,pRHS22 },
            });

            Matrix<double> pmatrix = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {0,0,1 },
                {0,1,0 },
                {1,0,0 },
            });

            Matrix<double> A = (pLHS*pmatrix).Solve(pRHS);
            Matrix<double> paramMatrix = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {2,3,-2,-3,0,0 },
                {2,3,0,0,-2,-3 },
                {0,0,2,3,-2,-3 },
                {1,1,0,0,0,0 },
                {0,0,1,1,0,0 },
                {0,0,0,0,1,1 },
            });
            Vector<double> b = Vector<double>.Build.Dense(new double[] { 0, 0, 0, 0, 0, 0 });
            Vector<double> qdotminus = Vector<double>.Build.Dense(new double[] { 0, 0, 0 });
            Vector<double> param = Vector<double>.Build.Dense(new double[] { 0, 0, 0, 0, 0, 0 });
            double anglepsi = (15.70*Math.PI)/32;
            //params
            double p0 = startq1;
            double p4 = startq2;
            double p8 = startq3;
            double p1;
            double p2;
            double p3;
            double p5;
            double p6;
            double p7;
            double p9;
            double p10;
            double p11;

            gait.vhc.phi1Parameters.Add("P0", p0);
            gait.vhc.phi1Parameters.Add("P1", 0);
            gait.vhc.phi1Parameters.Add("P2", 0);
            gait.vhc.phi1Parameters.Add("P3", 0);

            gait.vhc.phi2Parameters.Add("P4", p4);
            gait.vhc.phi2Parameters.Add("P5", 0);
            gait.vhc.phi2Parameters.Add("P6", 0);
            gait.vhc.phi2Parameters.Add("P7", 0);

            gait.vhc.phi3Parameters.Add("P8", p8);
            gait.vhc.phi3Parameters.Add("P9", 0);
            gait.vhc.phi3Parameters.Add("P10", 0);
            gait.vhc.phi3Parameters.Add("P11", 0);
            while (anglepsi > 0)
            {
                double a = 0.05;
                while(a < 1.5)
                {
                    double enddq1 = -a * Math.Cos(anglepsi);
                    double enddq3 = -a * Math.Sin(anglepsi);
                    
                    double enddq2 = -0.2;
                    while(enddq2< 0.2)
                    {
                        qdotminus[0] = enddq1;
                        qdotminus[1] = enddq2;
                        qdotminus[2] = enddq3;
                        Vector<double> qdotplus = A * qdotminus;
                        paramMatrix[0, 2] = -2 * (enddq1 / enddq2);
                        paramMatrix[0, 3] = -3 * (enddq1 / enddq2);
                        paramMatrix[1, 0] = 2 * (enddq3 / enddq1);
                        paramMatrix[1, 1] = 3 * (enddq3 / enddq1);
                        paramMatrix[2, 4] = -2 * (enddq2 / enddq3);
                        paramMatrix[2, 5] = -3 * (enddq2 / enddq3);
                        p1 = -0.05;
                        while(p1 > -2)
                        {
                            p5 = qdotplus[1] / qdotplus[0] * p1;
                            p9 = qdotplus[2] / qdotplus[0] * p1;

                            b[0] = -p1 + (enddq1 / enddq2) * p5;
                            b[1] = -p1* (enddq3 / enddq1) + p9;
                            b[2] = -p5 + (enddq2 / enddq3) * p9;
                            b[3] = -2 * p0 - p1;
                            b[4] = -p5;
                            b[5] = -2 * p8 - p9;
                            param = paramMatrix.Solve(b);
                            //param = paramMatrix.Inverse() * b;
                            gait.vhc.phi1Parameters["P1"] = p1;
                            gait.vhc.phi1Parameters["P2"] = param[0];
                            gait.vhc.phi1Parameters["P3"] = param[1];

                            gait.vhc.phi2Parameters["P5"] = p5;
                            gait.vhc.phi2Parameters["P6"] = param[2];
                            gait.vhc.phi2Parameters["P7"] = param[3];

                            gait.vhc.phi3Parameters["P9"] = p9;
                            gait.vhc.phi3Parameters["P10"] = param[4];
                            gait.vhc.phi3Parameters["P11"] = param[5];

                            bool neg = true;
                            for(int i = 0; i < 40; i++)
                            {
                                double val = gait.vhc.evalAlpha(0.025 * i);
                                Console.WriteLine("");
                                if (val > 1)
                                {
                                    neg = false;
                                    break;
                                }
                            }
                            if (neg)
                            {
                                using (StreamWriter file =
                                    new System.IO.StreamWriter(@"../../../foundgaits.txt", true))
                                {
                                    file.WriteLine(anglepsi);
                                    file.WriteLine(enddq2);
                                    file.WriteLine(p1);
                                    file.WriteLine(a);
                                    file.WriteLine(p0);
                                    file.WriteLine(p1);
                                    file.WriteLine(param[0]);
                                    file.WriteLine(param[1]);
                                    file.WriteLine(p4);
                                    file.WriteLine(p5);
                                    file.WriteLine(param[2]);
                                    file.WriteLine(param[3]);
                                    file.WriteLine(p8);
                                    file.WriteLine(p9);
                                    file.WriteLine(param[4]);
                                    file.WriteLine(param[5]);

                                    Console.WriteLine(anglepsi);
                                    Console.WriteLine(enddq2);
                                    Console.WriteLine(p1);
                                    Console.WriteLine(a);
                                    Console.WriteLine(p0);
                                    Console.WriteLine(p1);
                                    Console.WriteLine(param[0]);
                                    Console.WriteLine(param[1]);
                                    Console.WriteLine(p4);
                                    Console.WriteLine(p5);
                                    Console.WriteLine(param[2]);
                                    Console.WriteLine(param[3]);
                                    Console.WriteLine(p8);
                                    Console.WriteLine(p9);
                                    Console.WriteLine(param[4]);
                                    Console.WriteLine(param[5]);
                                    Console.WriteLine("");

                                }
                            }

                            p1 -= 0.05;
                        }

                        if(enddq2 == -0.01)
                        {
                            enddq2 += 0.02;
                        }
                        else
                        {
                            enddq2 += 0.01;
                        }
                        
                    }

                    a += 0.05;
                }


                anglepsi -= 0.01;
            }
        }
    }

    
}
