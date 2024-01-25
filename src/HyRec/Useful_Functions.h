
int Sign(double x)
{
  // Find sign of a number
  int result;
  if (x > 0.0)
  {
    result = 1;
  }
  else if (x < 0.0)
  {
    result = -1;
  }
  else
  {
    // maybe x is too close to 0?
    result = 0;
  }
  return result;
}

/*
int Same_Sign(double x, double y)
{
  int sx, sy, r;
  sx = Sign(x);
  sy = Sign(y);
  if (sx == sy)
  {
    r = 1;
  }
  else
  {
    r = 0;
  }
  return r;
}
*/

int Find_Index_1D(double *Tab, double x, int nx)
{
  // Find closest left element index
  // return -1 if not in range
  double x1, x2, x3;
  int id1, id2, id3, Stop, s1, s2, s3, idx, count;
  Stop = 0;
  id1 = 0;
  id3 = nx - 1;

  count = 0;
  while (Stop == 0)
  {
    count = count + 1;
    id2 = (int)round((((double)id1) + ((double)id3)) / 2.0);
    if (id3 == id1 + 1)
    {
      idx = id1;
      // printf("Stopping, id1 = %d, id2 = %d, id3 = %d\n", id1, id2, id3);
      Stop = 1;
    }

    x1 = Tab[id1];
    x2 = Tab[id2];
    x3 = Tab[id3];

    s1 = Sign(x - x1);
    s2 = Sign(x - x2);
    s3 = Sign(x - x3);
    if (s1 != s2)
    {
      id3 = id2;
    }
    else if (s2 != s3)
    {
      id1 = id2;
    }
    if ((s1 == s3) || (count > 50))
    {
      idx = -1; // this means x is not in range
      Stop = 1;
    }
  }
  return idx;
}

double Interp_1D(double x, double *x_axis, double *y_axis, int nx, int Overflow_Handle, int Interp_Method)
{
  /* Find value of y at x
  Overflow_Handle : what to do if x is not in x_axis
                    0 : raise error and exit
                    1 : give nearest value
                    2 : give 0
  Interp_Method : Interpolation method
                  0 : linear
                  1 : linear in log space
  */
  int id1, id2;
  double x1, x2, y1, y2, x_, r;
  id1 = Find_Index_1D(x_axis, x, nx);
  if (id1 == -1)
  {
    if (Overflow_Handle == 1)
    {
      if (x < x_axis[0])
      {
        r = y_axis[0];
      }
      else
      {
        r = y_axis[nx - 1];
      }
    }
    else if (Overflow_Handle == 2)
    {
      r = 0;
    }
    else
    {
      printf("Error from Interp_1D: x is not in range.\n");
      exit(1);
    }
  }
  else
  {
    id2 = id1 + 1;
    if (Interp_Method == 0)
    {
      x1 = x_axis[id1];
      x2 = x_axis[id2];
    }
    else if (Interp_Method == 1)
    {
      x1 = log(x_axis[id1]);
      x2 = log(x_axis[id2]);
      x_ = log(x);
    }
    else
    {
      printf("Error from Interp_1D: unkown interpolation method\n");
      exit(1);
    }

    y1 = y_axis[id1];
    y2 = y_axis[id2];
    r = (y2 - y1) / (x2 - x1) * (x_ - x1) + y1;

    //printf("x_ = %f, x1 = %f, x2 = %f, y1 = %f, y2 = %f\n", x_, x1, x2, y1, y2);
  }

  return r;
}

double Read_Tab_2D(int yid, int xid, int ny, int nx, double *Tab)
{
  // A way to read 2D table, I can't get pointer to work
  int id;
  if ((xid > nx - 1) || (yid > ny - 1))
  {
    printf("Error: array index overflow.\n");
    exit(1);
  }
  id = yid * nx + xid;
  return Tab[id];
}

double Interp_2D(double m, double z, double *m_axis, double *z_axis, double *Tab, int nm, int nz, int Overflow_Handle)
{
  /* Interpolate (in log) from a 2D table
   -- inputs --
   m : m target
   z : z target
   m_axis : m axis pointer
   z_axis : z axis pointer
   Tab : Data Table, index : Tab[z_id][m_id]
   nm : m axis size
   nz : z axis size
   Overflow_Handle : decide what to do when m or z is not in range
                     0 : Raise error and exit
                     1 : return 0 when z is not in range, if m is also not in range then raise error
  */
  int mid1, mid2, zid1, zid2;
  double lm, lm1, lm2, lz, lz1, lz2, f11, f12, f21, f22, f1, f2, F1, F2, f;
  mid1 = Find_Index_1D(m_axis, m, nm);
  zid1 = Find_Index_1D(z_axis, z, nz);
  // printf("debug: m = %E,  z = %E\n", m, z);

  if ((mid1 == -1) || (zid1 == -1))
  {
    // Overflow detected
    if (Overflow_Handle == 0)
    {
      printf("Error: m or z not in range, exitting, debugging info:.\n");
      printf("m1 = %E, m = %E, m2 = %E\n", m_axis[mid1], m, m_axis[mid1 + 1]);
      printf("z1 = %E, z = %E, z2 = %E\n", z_axis[zid1], z, z_axis[zid1 + 1]);
      printf("min(m) = %E, m = %E, max(m) = %E\n", m_axis[0], m, m_axis[nm - 1]);
      printf("min(z) = %E, z = %E, max(z) = %E\n", z_axis[0], z, z_axis[nz - 1]);
      exit(1);
    }
    else if ((zid1 == -1) && (Overflow_Handle == 1))
    {
      if (mid1 == -1)
      {
        printf("Error: m not in range, exitting, debugging info:.\n");
        printf("m1 = %E, m = %E, m2 = %E\n", m_axis[mid1], m, m_axis[mid1 + 1]);
        printf("z1 = %E, z = %E, z2 = %E\n", z_axis[zid1], z, z_axis[zid1 + 1]);
        printf("min(m) = %E, m = %E, max(m) = %E\n", m_axis[0], m, m_axis[nm - 1]);
        printf("min(z) = %E, z = %E, max(z) = %E\n", z_axis[0], z, z_axis[nz - 1]);
        exit(1);
      }
      else
      {
        // Our z axis is bounded
        return 0.0;
      }
    }
  }

  mid2 = mid1 + 1;
  zid2 = zid1 + 1;

  lm = log10(m);
  lm1 = log10(m_axis[mid1]);
  lm2 = log10(m_axis[mid2]);

  lz = log10(z);
  lz1 = log10(z_axis[zid1]);
  lz2 = log10(z_axis[zid2]);

  f11 = Read_Tab_2D(zid1, mid1, nz, nm, Tab);
  f12 = Read_Tab_2D(zid1, mid2, nz, nm, Tab);
  f21 = Read_Tab_2D(zid2, mid1, nz, nm, Tab);
  f22 = Read_Tab_2D(zid2, mid2, nz, nm, Tab);

  // fix m1
  f1 = f11;
  f2 = f21;
  F1 = (f2 - f1) / (lz2 - lz1) * (lz - lz1) + f1;
  // fix m2
  f1 = f12;
  f2 = f22;
  F2 = (f2 - f1) / (lz2 - lz1) * (lz - lz1) + f1;

  f = (F2 - F1) / (lm2 - lm1) * (lm - lm1) + F1;

  return f;
}

void Fill_Linear_Array(double xmin, double xmax, double *array, int nx, int method)
{
  /* Fill array with nx element (log) linearly spaced between [xmin, xmax]
  -- inputs --
  method : filling method
           0 : linear
           1 : linear log
   */
  double dx, x, x1, x2;
  int id;
  if (method == 0)
  {
    x1 = xmin;
    x2 = xmax;
    dx = (x2 - x1) / ((double)(nx - 1));
  }
  else if (method == 1)
  {
    x1 = log(xmin);
    x2 = log(xmax);
    dx = (x2 - x1) / ((double)(nx - 1));
  }
  x = x1;
  for (id = 0; id < nx; id++)
  {
    array[id] = exp(x);
    x += dx;
  }
}