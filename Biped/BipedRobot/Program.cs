using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BipedRobot{
    static class Program{
        static void Main(string[] args)
        {
            Biped biped = new Biped (@"../../basic_randomized_set.xml");
            integration.run (ref biped);
            plotting.plotStates(biped);
            gaitSearch.run(ref biped);

        }
    }
}
