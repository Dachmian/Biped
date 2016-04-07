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
            IntegrationFullDynamics.run(ref biped);
            plotting.plotStates(biped);
            gaitSearch.run(ref biped);
            //test();
            
        }

        static void test()
        {
            Expression test = Infix.ParseOrUndefined("a + b*x");
            Expression test2 = Infix.ParseOrUndefined("r + p*y^2");
            Expression result = Structure.Substitute("x", test2, test);
            Console.WriteLine(Infix.Format(result));
        }
    }
}
