using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace arma
{
    class Program
    {
        static void Main(string[] args)
        {
            arma test = new arma();
            test.dataInput("data_wind.txt");
            /* test.Diff();
             test.modelSpecification();
             test.foreCast(16);
             test.revDiff();*/
            test.armaProgress();
        }
    }
}
