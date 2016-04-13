using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Symbolics;

namespace BipedRobot{
    static class Program{
        static void Main(string[] args)
        {
            Biped biped = new Biped (@"../../basic_randomized_set.xml");
            //IntegrationFullDynamics.run(ref biped);
            //plotting.plotStates(biped);
            gaitSearch.run(ref biped);
            //test2(biped);
            
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
            
            int numberOfPoints = 8;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);
            double[] p = new double[30] { 1.002677806583060, 0.973626714059529, 0.974818663249441, 0.975161177554158, 0.973966193029470, 0.969903727868665, 0.960170394076285, 0.937599240792406, 0.888098904754680, 0.911104026678500, 0.986576776557014, 0.983218078674456, 0.987832792775939, 0.992391672194564, 0.996558129289077, 0.999730534966211, 1.000582664493830, 0.995792455591024, 0.975884520996868, 0.997017519177504, 0.898907034694419, 0.885021091002188, 0.891598071941415, 0.898378848358458, 0.905364451193932, 0.912687323064369, 0.920607258611730, 0.929670828893752, 1.100905857568060, 1.111317288007510 };
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
            Console.WriteLine(gait.impactFirstLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactSecondLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));
            Console.WriteLine(gait.impactThirdLine(gait.gaitParam.theta0, gait.gaitParam.thetaT));

            gaitSearch.setAndVerifyParameters(ref gait);
            Console.WriteLine(gait.gaitParam.dtheta0);
            Console.WriteLine(gait.gaitParam.dthetaT);
            
        }
    }
}
