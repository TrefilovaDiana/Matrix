using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;

namespace Matrix
{
    class Program
    {
        private class Pair //Класс пар для задания минимальных и максимальных значений элементов матриц
        {
            public int Max { get; } //Максимальное значение
            public int Avg { get; } //"Выравнивание" (если выравнивание 1000, то при генерации после получения r.Next(max) будет вычитаться 1000)
            public Pair(int avg, int max)
            {
                Max = max;
                Avg = avg;
            }
        }

        static void Read_Matrix(out int[,] A, string F_name) //Считывание матрицы из файла
        {
            using (StreamReader sr = new StreamReader(F_name))
            {
                string s = sr.ReadLine();
                int n = int.Parse(s);
                A = new int[n, n];
                for (int i = 0; i < n; i++)
                {
                    s = sr.ReadLine();
                    for (int j = 0; j < n; j++)
                        A[i, j] = int.Parse(s.Split()[j]);
                }
            }
        }

        public static int[,] Sum_m(int[,] a, int[,] b) //Сумма матриц
        {
            int n = (int)Math.Sqrt(a.Length);
            int[,] c = new int[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    c[i, j] = a[i, j] + b[i, j];
            return c;
        }

        public static int[,] Dif_m(int[,] a, int[,] b) //Сумма матриц
        {
            int n = (int)Math.Sqrt(a.Length);
            int[,] c = new int[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    c[i, j] = a[i, j] - b[i, j];
            return c;
        }

        public static int[,] Get_part(int i_from, int i_to, int j_from, int j_to, int[,] a) //Копирование части матрицы
        {
            int n = (int)Math.Sqrt(a.Length);
            int[,] c = new int[n / 2, n / 2];
            for (int i = i_from; i < i_to; i++)
                for (int j = j_from; j < j_to; j++)
                    c[i - i_from, j - j_from] = a[i, j];
            return c;
        }

        public static void Set_part(int i_from, int i_to, int j_from, int j_to, int[,] a, int[,] b) //Сохраняем часть матрицы в итоговую матрицу
        {
            for (int i = i_from; i < i_to; i++)
                for (int j = j_from; j < j_to; j++)
                    a[i, j] = b[i - i_from, j - j_from];
        }

        static void Generate(out int[,] A, int n, Pair range) //Генерация матриц
        {
            Random r = new Random();
            A = new int[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    A[i, j] = r.Next(range.Max) - range.Avg;

        }

        static void Generate_strassen(out int[,] A, int n, Pair range) //Генерация матриц для алгоритма Штрассена
        {
            int n_total = (int)Math.Pow(2, Math.Ceiling(Math.Log(n) / Math.Log(2))); //Берем Log(n)/Lon(2) округляем вверх и получаем ближайшую степень двойки
            Generate(out A, n_total, range); //Выполняем стандартную генерацию
        }

        static void Multiply(int[,] A, int[,] B, out int[,] C, int n) //Стандартное умножение матриц 
        {
            C = new int[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    C[i, j] = 0;
                    for (int k = 0; k < n; k++)
                        C[i, j] += A[i, k] * B[k, j];
                }
        }

        static void Strassen(int[,] A, int[,] B, out int[,] C) //Умножение с помощью алгоритма Штрассена
        {
            int n = (int)Math.Sqrt(A.Length);
            if (n == 1)
            {
                C = new int[n, n];
                C[0, 0] = A[0, 0] * B[0, 0];
                return;
            }

            //Разбиваем исходные матрицы на четверти
            int[,] a11 = Get_part(0, n / 2, 0, n / 2, A);
            int[,] a12 = Get_part(0, n / 2, n / 2, n, A);
            int[,] a21 = Get_part(n / 2, n, 0, n / 2, A);
            int[,] a22 = Get_part(n / 2, n, n / 2, n, A);

            int[,] b11 = Get_part(0, n / 2, 0, n / 2, B);
            int[,] b12 = Get_part(0, n / 2, n / 2, n, B);
            int[,] b21 = Get_part(n / 2, n, 0, n / 2, B);
            int[,] b22 = Get_part(n / 2, n, n / 2, n, B);


            //Далее - алгоритм Штрассена
            int[,] m1;
            Strassen(Sum_m(a11, a22), Sum_m(b11, b22), out m1);
            int[,] m2;
            Strassen(Sum_m(a21, a22), b11, out m2);
            int[,] m3;
            Strassen(a11, Dif_m(b12, b22), out m3);
            int[,] m4;
            Strassen(a22, Dif_m(b21, b11), out m4);
            int[,] m5;
            Strassen(Sum_m(a11, a12), b22, out m5);
            int[,] m6;
            Strassen(Dif_m(a21, a11), Sum_m(b11, b12), out m6);
            int[,] m7;
            Strassen(Dif_m(a12, a22), Sum_m(b21, b22), out m7);

            C = new int[n, n];
            Set_part(0, n / 2, 0, n / 2, C, Sum_m(Dif_m(Sum_m(m1, m4), m5), m7));
            Set_part(0, n / 2, n / 2, n, C, Sum_m(m3, m5));
            Set_part(n / 2, n, 0, n / 2, C, Sum_m(m2, m4));
            Set_part(n / 2, n, n / 2, n, C, Sum_m(Sum_m(Dif_m(m1, m2), m3), m6));
        }

        static void Print_m(int[,] A, int n) //Вывод матрицы (использовался только для отладки и тестов)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    Console.Write(A[i, j] + " ");
                Console.WriteLine();
            }
        }

        //static List<int> Generate_Size(int Max_size) //Заполняем список размерностей
        //{
            //List<int> A = new List<int>();
            //A.Add(1);
            //int i = 32;
            //while (i <= Max_size)
            //{
                //A.Add(i);
                //i += 32;
            //}
            //return A;
        //}

        static double Check_Average_Time()
        {
            int i = 0, count = 100000000;
            Stopwatch sw = new Stopwatch();
            sw.Start();
            while (i < count) //count сравнений + count инкрементов
            {
                i++;
            }
            sw.Stop();
            return (double)sw.ElapsedMilliseconds / count;
        }

        static void Main(string[] args)
        {
            Stopwatch sw = new Stopwatch();
            StreamWriter streamWriter = new StreamWriter("output1.txt");
            int Max_size = 512;
            int[,] A, B, C;
            int count_test = 100; //Количество тестов
            double time = 0;
            double time_total = 0;

            List<int> I = new List<int>();
            I.Add(1);
            I.Add(2);
            I.Add(4);
            I.Add(8);
            I.Add(16);
            I.Add(32);

            //Generate_Size(Max_size);

            List<Pair> Data = new List<Pair>();
            Data.Add(new Pair(0, 1000));
            Data.Add(new Pair(-1000, 2000));
            Data.Add(new Pair(-1000, 1000));

            double avg = Check_Average_Time();

            streamWriter.WriteLine("Среднее время выполнения условной операции: " + avg + " мс.");
            streamWriter.WriteLine();
            foreach (Pair range in Data)
            {
                streamWriter.WriteLine("Значения элементов матриц от " + range.Avg + " до " + (range.Max + range.Avg));
                streamWriter.WriteLine();
                streamWriter.WriteLine("Стандартное умножение матриц:");
                streamWriter.WriteLine();
                streamWriter.WriteLine("{0:15}|{1:16}", "Размер матрицы", "Время выполнения");
                foreach (int item in I) //Выполняем count_test тестов для каждой размерности в I
                {
                    time = 0;
                    for (int i = 0; i < count_test; i++)
                    {
                        //Генерируем матрицы
                        Generate(out A, item, range);
                        Generate(out B, item, range);

                        //Сбрасываем таймер
                        sw.Restart();
                        //Выполняем умножение
                        Multiply(A, B, out C, item);
                        //Останавливаем таймер
                        sw.Stop();
                        //Суммируем затраченное время
                        time += sw.ElapsedMilliseconds;
                    }
                    //Находим среднее время выполнения теста
                    time_total = time / count_test;
                    //Записываем результат в файл
                    streamWriter.WriteLine("{0,14} {1,16}", item.ToString(), time_total);
                }

                streamWriter.WriteLine();


                streamWriter.WriteLine("Алгоритм Штрассена:");
                streamWriter.WriteLine("{0:15}|{1:16}", "Размер матрицы", "Время выполнения");
                foreach (int item in I) //Выполняем count_test тестов для каждой размерности в I
                {
                    time = 0;
                    for (int i = 0; i < count_test; i++)
                    {
                        //Генерируем матрицы
                        Generate_strassen(out A, item, range);
                        Generate_strassen(out B, item, range);

                        //Сбрасываем таймер
                        sw.Restart();
                        //Выполняем умножение
                        Strassen(A, B, out C);
                        //Останавливаем таймер
                        sw.Stop();
                        //Суммируем затраченное время
                        time += sw.ElapsedMilliseconds;
                    }
                    //Находим среднее время выполнения теста
                    time_total = time / count_test;
                    //Записываем результат в файл
                    streamWriter.WriteLine("{0,14} {1,16}", item.ToString(), time_total);
                }
                streamWriter.WriteLine();
            }
            streamWriter.Close();
        }
    }
}
