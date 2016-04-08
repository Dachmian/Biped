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
            
            int numberOfPoints = 5;
            BRgait gait = new BRgait(biped.param, numberOfPoints);
            BezierCurve brCrv = new BezierCurve(numberOfPoints, gait);
            
            Tuple<Dictionary<string, FloatingPoint>, string> tuple = brCrv.phi1ToString();
            gait.vhc.phi1 = Infix.ParseOrUndefined(tuple.Item2);
            tuple.Item1["P0"] = -0.397706480340997;
            tuple.Item1["P1"] = -1.832099220573480;
            tuple.Item1["P2"] = -1.772153266286530;
            tuple.Item1["P3"] = 0.131396187920137;
            tuple.Item1["P4"] = 1.935281527481110;
            gait.vhc.phi1Parameters = tuple.Item1;

            tuple = brCrv.phi2ToString();
            gait.vhc.phi2 = Infix.ParseOrUndefined(tuple.Item2);
            tuple.Item1["P5"] = -1.007680385447250;
            tuple.Item1["P6"] = 0.166771260567510;
            tuple.Item1["P7"] = 4.708400907030840;
            tuple.Item1["P8"] = -3.867479905112620;
            tuple.Item1["P9"] = 0.940408914676057;
            gait.vhc.phi2Parameters = tuple.Item1;
            //{0.348867076182787,1.734452411850560,1.745040104384210,0.323596856033392,-2.075955418952420}
            tuple = brCrv.phi3ToString();
            gait.vhc.phi3 = Infix.ParseOrUndefined(tuple.Item2);
            tuple.Item1["P10"] = 0.348867076182787;
            tuple.Item1["P11"] = 1.734452411850560;
            tuple.Item1["P12"] = 1.745040104384210;
            tuple.Item1["P13"] = 0.323596856033392;
            tuple.Item1["P14"] = -2.075955418952420;
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
        }
    }
}
