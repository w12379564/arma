using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;


namespace arma
{
    class arma
    {
        int pOrder;
        int qOrder;
        int p;
        int q;
        double[] coefficientAR;
        double[] coefficientMA;
        double[] data;
        double[] eData;//原序列残差
        double[] dataNew;
        double[] foreValue;
       
        double[] dataCopy;

        public double[] outputData;

        Random ranGauss = new Random(100);

        /*public arma(int arOrder,int maOrder)
        {
            p = arOrder;
            q = maOrder;
            coeAR = new double[p];
            coeMA = new double[q];

        }*/

        public void dataInput(String path)
        {
            ArrayList dataList = new ArrayList();
            StreamReader sr = new StreamReader(path);

            String line;
            while ((line = sr.ReadLine()) != null)
            {
                dataList.Add(Convert.ToDouble(line));
            }

            data = new double[dataList.Count];

            for (int i = 0; i < dataList.Count; i++)
            {
                data[i] = (double)dataList[i];
            }


            dataCopy = new double[data.Length];
            data.CopyTo(dataCopy, 0);

            Console.WriteLine("Data loaded.");

        }

        public void dataInput(double[] dataList)
        {

            data = new double[dataList.Length];

            for (int i = 0; i < dataList.Length; i++)
            {
                data[i] = dataList[i];
            }


            dataCopy = new double[data.Length];
            data.CopyTo(dataCopy, 0);

            Console.WriteLine("Data loaded.");
        }

        public void armaProgress()
        {
            Diff();
            modelSpecification();
            foreCast(16);
            revDiff();

            outputData = new double[foreValue.Length];
            for(int i=0;i<foreValue.Length;i++)
            {
                outputData[i] = foreValue[i];
            }

        }

        public double residualFunc(double[] coeRow)
        {
            double[] e = new double[data.Length];
            double sumSquare = 0;
            int arIndex;
            int maIndex;

            double[] coeAR = new double[p];
            double[] coeMA = new double[q];

            for (int i = 0; i < p; i++)
            {
                coeAR[i] = coeRow[i];
            }

            for (int i = 0; i < q; i++)
            {
                coeMA[i] = coeRow[i + p];
            }

            for (int t = 0; t < data.Length; t++)
            {
                arIndex = t - 1;
                maIndex = t - 1;

                e[t] = e[t] + data[t];

                while (arIndex >= 0 && (t - arIndex) <= p)
                {
                    e[t] = e[t] + coeAR[t - arIndex - 1] * data[arIndex];
                    arIndex -= 1;
                }

                while (maIndex >= 0 && (t - maIndex) <= q)
                {
                    e[t] = e[t] - coeMA[t - maIndex - 1] * e[maIndex];
                    maIndex -= 1;
                }

            }

            for (int i = 0; i < e.Length; i++)
            {
                sumSquare += e[i] * e[i];
            }

            return sumSquare;
        }

        double[] getResidual(int arOrder, int maOrder, double[] coeAR, double[] coeMA)
        {
            int arIndex;
            int maIndex;
            double[] e = new double[data.Length];

            for (int t = 0; t < data.Length; t++)
            {
                arIndex = t - 1;
                maIndex = t - 1;

                e[t] = e[t] + data[t];

                while (arIndex >= 0 && (t - arIndex) <= arOrder)
                {
                    e[t] = e[t] + coeAR[t - arIndex - 1] * data[arIndex];
                    arIndex -= 1;
                }

                while (maIndex >= 0 && (t - maIndex) <= maOrder)
                {
                    e[t] = e[t] - coeMA[t - maIndex - 1] * e[maIndex];
                    maIndex -= 1;
                }

            }
            return e;
        }


        public void modelSpecification()
        {

            Console.WriteLine("Model Specification.");

            int maxp = 3;
            int maxq = 3;

            int[] pArray = new int[maxp * maxq];
            int[] qArray = new int[maxp * maxq];
            double[] AIC = new double[maxp * maxq];
            double[] solution;
            double[][] coeARlist = new double[maxp * maxq][];
            double[][] coeMAlist = new double[maxp * maxq][];
            double[] e;
            double varValue;
            int minAICIndex;

            int count = 0;

            for (p = 1; p <= maxp; p++)
            {

                for (q = 1; q <= maxq; q++)
                {

                    Console.Write("Testing ARMA({0},{1}) model......", p, q);

                    pArray[count] = p;
                    qArray[count] = q;
                    PSO pso_calc = new PSO(p + q);
                    pso_calc.objFunc = new PSO.objHandler(this.residualFunc);
                    solution = pso_calc.pso_iteration();
                    coeARlist[count] = new double[p];
                    coeMAlist[count] = new double[q];

                    for (int i = 0; i < p; i++)
                    {
                        coeARlist[count][i] = solution[i];
                    }

                    for (int i = 0; i < q; i++)
                    {
                        coeMAlist[count][i] = solution[p + i];
                    }

                    e = getResidual(p, q, coeARlist[count], coeMAlist[count]);

                    varValue = var(e);

                    AIC[count] = e.Length * Math.Log10(varValue) + 2 * (p + q + 2);

                    count++;

                    Console.WriteLine("Done");

                }
            }

            minAICIndex = minIndex(AIC);
            pOrder = pArray[minAICIndex];
            qOrder = qArray[minAICIndex];

            coefficientAR = coeARlist[minAICIndex];
            coefficientMA = coeMAlist[minAICIndex];

            Console.WriteLine("The model is ARMA({0},{1}) model.\n", pOrder, qOrder);
            Console.Write("The expression is:\na[t]");
            for (int i = 0; i < pOrder; i++)
            {
                Console.Write(" + {0} * a[t-{1}]", coefficientAR[i], i + 1);
            }
            Console.Write(" = e[t]");
            for (int i = 0; i < qOrder; i++)
            {
                Console.Write(" + {0} * e[t-{1}]", coefficientMA[i], i + 1);
            }
            Console.WriteLine("");

        }

        double armaExpression(double[] e, int t)
        {
            double at = 0;

            int arIndex = t - 1;
            int maIndex = t - 1;

            while (arIndex >= 0 && (t - arIndex) <= pOrder)
            {
                at -= coefficientAR[t - arIndex - 1] * dataNew[arIndex];
                arIndex--;
            }

            at += e[t];

            while (maIndex >= 0 && (t - maIndex) <= qOrder)
            {
                at += coefficientMA[t - maIndex - 1] * e[maIndex];
                maIndex--;
            }

            return at;

        }

        public double[] foreCast(int num)
        {

            Console.WriteLine("\nThe next {0} forecast value is:", num);

            Random ran = new Random();
            double[] eFore = new double[num];
            double[] eNew;

            dataNew = new double[data.Length + num];
            foreValue = new double[num];

            for (int i = 0; i < dataNew.Length; i++)
            {
                if (i < data.Length)
                {
                    dataNew[i] = data[i];
                }
                else
                {
                    dataNew[i] = 0;
                }
            }

            eData = getResidual(pOrder, qOrder, coefficientAR, coefficientMA);

            double mean_e = mean(eData);
            double var_e = var(eData);

            for (int i = 0; i < num; i++)//需要用白噪声程序替换
            {
                eFore[i] = Gaussian(mean_e, Math.Sqrt(var_e));
            }

            eNew = new double[eData.Length + eFore.Length];

            for (int i = 0; i < eNew.Length; i++)
            {
                if (i < eData.Length)
                {
                    eNew[i] = eData[i];
                }
                else
                {
                    eNew[i] = eFore[i - eData.Length];
                }
            }

            for (int t = 0; t < foreValue.Length; t++)
            {
                foreValue[t] = armaExpression(eNew, t + data.Length);
                dataNew[t + data.Length] = foreValue[t];
                //Console.WriteLine("{0}", foreValue[t]);
            }



            return foreValue;
        }

        double mean(double[] x)
        {
            double sum = 0;

            for (int i = 0; i < x.Length; i++)
            {
                sum += x[i];
            }

            return sum / x.Length;
        }

        double var(double[] x)
        {
            double meanValue = mean(x);
            double sum = 0;

            for (int i = 0; i < x.Length; i++)
            {
                sum += (x[i] - meanValue) * (x[i] - meanValue);
            }

            return sum / x.Length;
        }

        int minIndex(double[] x)
        {
            for (int i = 0; i < x.Length; i++)
            {
                if (x[i] == x.Min())
                {
                    return i;
                }
            }
            return -1;
        }


        double Gaussian(double mu, double sigma)
        {
            
                  
            double u1 = ranGauss.NextDouble();
            double u2 = ranGauss.NextDouble();

            double z = Math.Sqrt(-2 * Math.Log(u1)) * Math.Cos(2 * Math.PI * u2);

            return mu + sigma * z;

        }
        
        public void Diff()
        {
            data = new double[data.Length - 1];
            for (int i = 0; i < data.Length; i++)
            {
                data[i] = dataCopy[i + 1] - dataCopy[i];
            }

        }

        public void revDiff()
        {

            foreValue[0] += dataCopy[dataCopy.Length - 1];

            for (int i = 1; i < foreValue.Length; i++)
            {
                foreValue[i] = foreValue[i - 1] + foreValue[i];
                Console.WriteLine("{0}", foreValue[i]);
            }
        }


    }
}
