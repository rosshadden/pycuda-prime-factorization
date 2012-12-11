#ifndef __qs_h
#define __qs_h

typedef struct factor_base_element
{
  int p, logp, y, z;
  int sqrtn1, sqrtn2;
} factor_base_element;

inline void *Calloc(int elements, int length)
{
  void *p;

  if (!(p = calloc(elements, length))) {
     fprintf(stderr,"Memory allocation failed for %d elements of length %d\n", elements, length);
     fflush(stderr);
     abort();
     return NULL;
   }
   return p;
}

/* Use the Shanks-Tonelli algorithm to solve x^2 = a (mod q) */
inline int mpz_sqrtm(mpz_ptr rop, mpz_t a, mpz_t q)
{
  mpz_t g, temp, t, gInv, qDiv, h, b;
  int i, s, e, y;

  if (mpz_legendre(a, q) == -1)
    return 0;

  mpz_init(g); mpz_init(temp);

  mpz_sub_ui(temp, q, 1);

  while (mpz_legendre( g, q ) != -1)
    {
      mpz_random(g, 2);
      mpz_mod(g, g, temp);
      mpz_add_ui(g, g, 1);
    }

  mpz_init_set(t, q);
  mpz_sub_ui(t, t, 1);
  s = mpz_scan1(t, 0);
  mpz_tdiv_q_2exp(t, t, s);

  e = 0;

  mpz_init(gInv);
  if (!mpz_invert(gInv, g, q))
    return 0;

  mpz_init(qDiv);
  mpz_init(h);
  for (i = 2; i <= s; i++)
    {
      mpz_powm_ui(temp, gInv, e, q);
      mpz_mul(h, a, temp);
      mpz_sub_ui(temp, q, 1);
      mpz_tdiv_q_2exp(qDiv, temp, i);
      mpz_powm(temp, h, qDiv, q);
      if (mpz_cmp_ui(temp, 1 ) != 0)
	{
	  y = 1 << (i - 1);
	  e += y;
	}
    }

  mpz_powm_ui(temp, gInv, e, q);
  mpz_mul(h, a, temp);

  mpz_init(b);
  mpz_add_ui(t, t, 1);
  mpz_tdiv_q_2exp(t, t, 1);
  mpz_powm(h, h, t, q);
  mpz_powm_ui(g, g, (int)e/2, q);
  mpz_mul(b, g, h);
  mpz_mod(b, b, q);
  mpz_set(rop, b);

  mpz_clear(g); mpz_clear(temp); mpz_clear(t); mpz_clear(gInv); mpz_clear(qDiv); mpz_clear(h); mpz_clear(b);
  return 1;
}

/* A wrapper for the shanks-tonelli algorithm */
int tonelli(int a, int p)
{
  int resi;
  mpz_t b, q, res;
  mpz_init(b); mpz_init(q); mpz_init(res);
  mpz_set_ui(b, a);
  mpz_set_ui(q, p);
  mpz_sqrtm(res, b, q);
  resi = mpz_get_ui(res);
  mpz_clear(b); mpz_clear(q); mpz_clear(res);
  return resi;
}

/* Augment a matrix in GF(2) with the identity matrix and perform gaussian elimination */
int gauss_eliminate(int r, mpz_t *a, int extra)
{
  int i, j, k, l, found;

  /* Initialize the identity matrix */
  for (i = 0; i < r; i++)
    {
      mpz_setbit(a[i], i + r - extra);
    }

  /* Iterate through all rows */
  for (i = 0; i < r; i++)
    {
      /* Does this row contain only zeros */
      if (mpz_cmp_ui(a[i], 0) == 0)
	{
	  /* If there is an all zero row below this one then swap the rows */
	  found = 0;
	  for (j = i+1; j < r; j++)
	    if (mpz_cmp_ui(a[j], 0) != 0)
	      {
		found = 1;
		break;
	      }
	  if (found == 0)
	    return 0;
	  else {
	    mpz_swap(a[i], a[j]);
	  }
	}

      /* Find the index of the first non-zero entry in this row and if there is a non-zero entry in the matrix below and to the left then swap */
      found = 1;
      while (found)
	{
	  found = 0;
	  j = mpz_scan1(a[i], 0);
	  if (j > 0)
	    {
	      for (k = i+1; k < r; k++)
		{
		  l = mpz_scan1(a[k], 0);
		  if (l < j)
		    {
		      mpz_swap(a[i], a[k]);
		      found = 1;
		      break;
		    }
		}
	    }
	}

      /* Change every other rows j'th entry to 0 */
      for (k = 0; k < r; k++)
	{
	  if (i != k)
	    {
	      l = mpz_scan1(a[k], 0);
	      if (l == j)
		{
		  mpz_xor(a[k], a[i], a[k]);
		}
	    }
	}
    }
  return 0;
}

#endif
