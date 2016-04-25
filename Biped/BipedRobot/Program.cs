using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Symbolics;
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
            findGaitsManually(biped);
            
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
            Console.WriteLine(Evaluate.Evaluate(_parameters, test).RealValue);
            
        }

        static void test2(Biped biped)
        {
            
            int numberOfPoints = 4;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);
            double[] p = { -1.12283194505, 4.79402021549, -0.761855899323, -4.30651425382, -2.43197378619, 2.44812379607, 2.64854178928, 0.613506850183, -4.84721295451, -0.909906754727, 3.50138134178, -4.46585242224 };
            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);

            for (int i = 0; i < p.Length; i++)
            {
                tuple.Item1["P" + i.ToString()] = p[i];

            }
            gait.vhc.phi1Parameters = tuple.Item1;
            
            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 0; i < p.Length; i++)
            {
                tuple.Item1["P" + i.ToString()] = p[i];

            }
            gait.vhc.phi2Parameters = tuple.Item1;
            
            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            for (int i = 0; i < p.Length; i++)
            {
                tuple.Item1["P" + i.ToString()] = p[i];

            }
            gait.vhc.phi3Parameters = tuple.Item1;

            gait.vhc.dphi1 = Infix.ParseOrUndefined(brCrv.dphi1ToString());
            gait.vhc.dphi2 = Infix.ParseOrUndefined(brCrv.dphi2ToString());
            gait.vhc.dphi3 = Infix.ParseOrUndefined(brCrv.dphi3ToString());

            gait.vhc.ddphi1 = Infix.ParseOrUndefined(brCrv.ddphi1ToString());
            gait.vhc.ddphi2 = Infix.ParseOrUndefined(brCrv.ddphi2ToString());
            gait.vhc.ddphi3 = Infix.ParseOrUndefined(brCrv.ddphi3ToString());

            Console.WriteLine(Infix.Format(gait.vhc.phi1));
            Console.WriteLine(Infix.Format(gait.vhc.phi2));
            Console.WriteLine(Infix.Format(gait.vhc.phi3));

            Console.WriteLine(Infix.Format(gait.vhc.dphi1));
            Console.WriteLine(Infix.Format(gait.vhc.dphi2));
            Console.WriteLine(Infix.Format(gait.vhc.dphi3));

            Console.WriteLine(Infix.Format(gait.vhc.ddphi1));
            Console.WriteLine(Infix.Format(gait.vhc.ddphi2));
            Console.WriteLine(Infix.Format(gait.vhc.ddphi3));

            Console.WriteLine(gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactSecondLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactThirdLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));

            //gaitSearch.setAndVerifyParameters(ref gait);
            //Console.WriteLine(gait.gaitParam.dtheta0);
            //Console.WriteLine(gait.gaitParam.dthetaT);
            
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




            //Console.WriteLine(gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            //Console.WriteLine(gait.impactSecondLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            //Console.WriteLine(gait.impactThirdLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));


            double[,] THETA = evalDthetaConstraint(gait, Math.Pow(dtheta0, 2), Math.Pow(dthetaT, 2));

            Phaseportrait plot = new Phaseportrait(THETA);
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
    }

    
}
