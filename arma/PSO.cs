using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace arma
{
    class PSO
    {
        private double minFitness = 9999999;
        private int minIndex = 1;
        private int n_var;
        private int n_particle;


        private double[,] particle;
        private double[,] velocity;
        private double[] v_max;

        private double[,] xLocalbest;
        private double[] fLocalbest;

        private double[] xGlobalbest;
        private double fGlobalbest;
        public delegate double objHandler(double[] x);
        //public delegate double objHandler(int p, int q, double[] coeAR, double[] coeMA);  //objetive function handler
        public objHandler objFunc;


        public PSO(int n)
        {
            n_var = n;
            n_particle = 20;
            particle = new double[n_particle, n_var];
            velocity = new double[n_particle, n_var];
            v_max = new double[n_var];

            xLocalbest = new double[n_particle, n_var];
            fLocalbest = new double[n_particle];

            xGlobalbest = new double[n_var];
            fGlobalbest = 9999999;

            //obj func = new obj();
            //objFunc = new objHandler(func.objFunc);

            //arma func = new arma();
            //objFunc = new objHandler(func.residualFunc);

            Random ran = new Random();

            for (int i = 0; i < n_var; i++)
            {
                v_max[i] = 0.1;  //ub of v
            }

            for (int i = 0; i < fLocalbest.Length; i++)
            {
                fLocalbest[i] = 9999999;
            }

            for (int i = 0; i < n_particle; i++)
            {
                for (int j = 0; j < n_var; j++)
                {
                    particle[i, j] = 10 * (ran.NextDouble() - 0.5);
                    velocity[i, j] = 0;
                }
            }

        }


        public double[] pso_iteration()
        {
            int t = 0;
            int t_last = 0;

            double[] score = new double[n_particle];
            double[] particle_row = new double[n_var];
            int[] betterIndex = new int[n_particle];
            double r1;
            double r2;
            Random ran = new Random();

            double w = 0.9;
            double c1 = 1;
            double c2 = 2;



            while (t - t_last < 100)
            {
                t = t + 1;

                for (int i = 0; i < n_particle; i++)
                {
                    for (int j = 0; j < n_var; j++)
                    {
                        particle_row[j] = particle[i, j];
                    }

                    score[i] = objFunc(particle_row);
                }

                betterIndex = find_better(score, fLocalbest);

                foreach (int index in betterIndex)
                {
                    if (index != -1)
                    {
                        fLocalbest[index] = score[index];
                        for (int j = 0; j < n_var; j++)
                        {
                            xLocalbest[index, j] = particle[index, j];
                        }
                    }
                }

                find_best(fLocalbest);

                if (t == 0 || minFitness + 0.01 < fGlobalbest)
                {
                    fGlobalbest = minFitness;

                    for (int i = 0; i < n_var; i++)
                    {
                        xGlobalbest[i] = particle[minIndex, i];
                    }
                }

                r1 = ran.NextDouble();
                r2 = ran.NextDouble();

                for (int i = 0; i < n_particle; i++)
                {
                    for (int j = 0; j < n_var; j++)
                    {
                        velocity[i, j] = w * velocity[i, j] + c1 * r1 * (xLocalbest[i, j] - particle[i, j]) + c2 * r2 * (xGlobalbest[j] - particle[i, j]);

                        if (velocity[i, j] > v_max[j])
                        {
                            velocity[i, j] = v_max[j];
                        }
                        if (velocity[i, j] < -v_max[j])
                        {
                            velocity[i, j] = -v_max[j];
                        }

                        particle[i, j] += velocity[i, j];
                    }
                }

                //Console.WriteLine("{0}\n", fGlobalbest);
            }

            return xGlobalbest;
        }


        void find_best(double[] x)
        {

            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] < minFitness)
                {
                    minFitness = x[i];
                    minIndex = i;

                }
            }
        }

        int[] find_better(double[] a, double[] b)
        {
            int[] betterIndex = new int[a.Length];
            for (int i = 0; i < a.Length; i++)
            {
                betterIndex[i] = -1;
            }

            for (int i = 0; i < a.Length; i++)
            {
                if (a[i] < b[i])
                {
                    betterIndex[i] = i;
                }
            }

            return betterIndex;
        }

        public void print()
        {
            Console.WriteLine("minimun x:{0} minimun f:{1}", xGlobalbest, fGlobalbest);
        }

        /*double objFunc(double[] x)
        {
            return Math.Pow(x[0], 2) - 4 * x[0] + 7;
        }*/

    }
}
