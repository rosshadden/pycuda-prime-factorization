#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <sys/time.h>
#include <unistd.h>
#include <time.h>
#include <signal.h>
#include <mpi.h>
#include "qs.h"

#define AFXTAG 1
#define AXBTAG 2
#define RESULTACKTAG 4
#define DTAG 8
#define NTAG 12
#define REQDTAG 16
#define REQNTAG 20
#define RESULTTAG 24
#define DIETAG 28
#define NLENGTHTAG 32

#define EXTRA 5
#define SIEVE_LENGTH 500000

#define MPI_BUFFER_SIZE 131072

/* Avoid a stupid warning by declaring this extern */
extern FILE * fmemopen (void *buf, size_t size, const char *opentype);

int threshold;
unsigned char *logfxs;
unsigned int factor_base_size;
factor_base_element *factor_base;
mpz_t *exponent_matrix, *exponent_vector, *axb, *afx;
mpz_t N, A, B, C, M, D;
double T;

/* Compute the sieve polynomial of x and store the result in y */
#define poly(x, y)         \
        mpz_mul(y, A, x);  \
        mpz_add(y, y, B);  \
        mpz_add(y, y, B);  \
        mpz_mul(y, y, x);  \
        mpz_add(y, y, C)

void send_mpz(int tag, int receiver, mpz_srcptr value, char *buffer)
{
  FILE *bp;
  int size;

  size = 4 + (mpz_size(value) * ((mp_bits_per_limb / 8) + (mp_bits_per_limb % 8 == 0 ? 0 : 1)));
  if (size > MPI_BUFFER_SIZE)
    {
      printf("error: MPI_BUFFER_SIZE (%d) is too small to hold %d bytes. Try to increase MPI_BUFFER_SIZE\n", MPI_BUFFER_SIZE, size);
    }

  bp = (FILE *)fmemopen(buffer, MPI_BUFFER_SIZE, "w");
  size = mpz_out_raw(bp, value);
  fclose(bp);

  MPI_Send(buffer, size, MPI_BYTE, receiver, tag, MPI_COMM_WORLD);
}

void buffer2mpz(mpz_ptr op, char *buffer)
{
  FILE *bp;

  bp = (FILE *)fmemopen(buffer, MPI_BUFFER_SIZE, "r");
  mpz_inp_raw(op, bp);
  fclose(bp);
}

/* Initialize the master node. x is the length of n */
int init_master(char *number, int x)
{
  int i;

  /* Return false if n has less than 10 digits */
  if (x < 9) return 0;

  /* Initialize numbers */
  mpz_init(N); mpz_init(M); mpz_init(D);
  mpz_set_str(N, number, 10);

  /* Set D to indicate that we have not found a suitable D yet */
  mpz_set_ui(D, 0);

  /* Compute the factor base size */
  if (x < 25)
    factor_base_size = 100;
  else
    factor_base_size = (int)(2.93 * (x * x) - 164.4 * x + 2455);

  x = 386 * (x * x) - 23209.3 * x + 352768;
  mpz_set_ui(M, x);

  /* Allocate space for the (ax + b) values */
  axb = (mpz_t *)Calloc(factor_base_size + EXTRA, sizeof(mpz_t));
  for (i = 0; i < factor_base_size + EXTRA; i++) mpz_init(axb[i]);

  /* Allocate space for the (a * f(x)) values */
  afx = (mpz_t *)Calloc(factor_base_size + EXTRA, sizeof(mpz_t));
  for (i = 0; i < factor_base_size + EXTRA; i++) mpz_init(afx[i]);

  /* Allocate space for the exponent matrix */
  exponent_matrix = (mpz_t *)Calloc(factor_base_size + EXTRA, sizeof(mpz_t));
  mpz_array_init((mpz_ptr)exponent_matrix, factor_base_size + EXTRA, factor_base_size * 2 + EXTRA);

  return 1;
}

/* We only call this function when N has been initialised. x is the length of n */
int init_slave(int x)
{
  /* Initialize the value to be factored */
  mpz_init(A); mpz_init(B); mpz_init(C); mpz_init(D); mpz_init(M);

  /* Calculate size of the factor base, length of the sieve interval and part of the threshold */
  if (x < 25)
    factor_base_size = 100;
  else
    factor_base_size = (int)(2.93 * (x * x) - 164.4 * x + 2455);
  T = 0.0268849 * x + 0.783929;
  x = 386 * (x * x) - 23209.3 * x + 352768;
  mpz_set_ui(M, x);

  /* Allocate space for the factor base */
  factor_base = (factor_base_element *)Calloc(factor_base_size, sizeof(factor_base_element));

  /* Allocate space for the logfx values */
  logfxs = (unsigned char *)Calloc(SIEVE_LENGTH, sizeof(unsigned char));

  /* Allocate space for the exponent vector */
  exponent_vector = (mpz_t *)Calloc(1, sizeof(mpz_t));
  mpz_array_init((mpz_ptr)exponent_vector, 1, factor_base_size);

  return 1;
}

void get_next_D()
{
  mpz_t tmp;

  mpz_init(tmp);

  if (mpz_cmp_ui(D, 0) == 0)
    {
      /* Set D to sqrt(sqrt(2N)/2) */
      mpz_mul_ui(D, N, 2);
      mpz_sqrt(D, D);
      mpz_cdiv_q(D, D, M);
      mpz_sqrt(D, D);
    }

  /* Find the next prime after D that is a quadratic residue modulo N such that D = 3 (mod 4) */
  mpz_nextprime(D, D);
  while ((mpz_legendre(N, D) != 1) || (mpz_tdiv_r_ui(tmp, D, 4) != 3))
    {
      mpz_nextprime(D, D);
    }

  mpz_clear(tmp);
}

/* Compute new values for A, B, C and store them in the global variables */
void generate_new_polynomial()
{
  mpz_t tmp, h0, h1, h2;

  mpz_init(tmp); mpz_init(h0); mpz_init(h1); mpz_init(h2);

  /* Get a new A */
  get_next_D();
  mpz_mul(A, D, D);

  /* Compute B */

  /* h0 = N^((D-3) / 4) */
  mpz_sub_ui(tmp, D, 3);
  mpz_cdiv_q_ui(tmp, tmp, 4);
  mpz_powm(h0, N, tmp, D);

  /* h1 = N^((D+1) / 4) */
  mpz_add_ui(tmp, D, 1);
  mpz_cdiv_q_ui(tmp, tmp, 4);
  mpz_powm(h1, N, tmp, D);

  /* h2 = (2*h1)^(-1) */
  mpz_mul_ui(tmp, h1, 2);
  mpz_invert(h2, tmp, D);

  /* tmp = (N - h1^2) / D */
  mpz_mul(tmp, h1, h1);
  mpz_sub(tmp, N, tmp);
  mpz_cdiv_q(tmp, tmp, D);

  /* h2 = (2*h1)^(-1) * [(N - h1^2) / D] */
  mpz_mul(h2, h2, tmp);

  /* B = h1 + h2*D mod A */
  mpz_mul(tmp, h2, D);
  mpz_add(B, h1, tmp);
  mpz_tdiv_r(B, B, A);

  /* Compute C */
  mpz_mul(C, B, B);
  mpz_sub(C, N, C);
  mpz_tdiv_q(C, C, A);
  mpz_neg(C, C);

  mpz_clear(tmp); mpz_clear(h0); mpz_clear(h1); mpz_clear(h2);
}

/* Solve f(x) = 0 (mod p) for each prime p in the factor base */
void initialize_polynomial()
{
  int i, j, c, sqroot, n;
  mpz_t tmp, inv, p;

  mpz_init(p); mpz_init(tmp); mpz_init(inv);

  sqroot = 0;
  for (i = 1; i < factor_base_size; i++)
    {
      mpz_set_ui(p, factor_base[i].p);
      n = mpz_tdiv_r_ui(tmp, N, factor_base[i].p);

      /* The two solutions to x^2 = N (mod p) - That is the two squareroots of N (mod p) */
      for (j = 0; j < 2; j++)
	{
	  sqroot = (j == 0 ? factor_base[i].sqrtn1 : factor_base[i].sqrtn2);
	  mpz_mul_ui(tmp, B, 2);
	  mpz_neg(tmp, tmp);
	  mpz_add_ui(tmp, tmp, 2 * sqroot);

	  /* Multiply by the inverse of (2A) (mod p) */
	  mpz_mul_ui(inv, A, 2);

	  /* Do it the trivial way if we can't find an inverse of (2A) (mod p) */
	  /* This is only the case if 2A = 0 (mod p) which will never occur for large values of A */
	  if (mpz_invert(inv, inv, p) == 0)
	    {
	      mpz_set_ui(inv, 0);
	      poly(inv, tmp);
	      mpz_tdiv_r_ui(tmp, tmp, factor_base[i].p);
	      while (mpz_cmp_ui(tmp, 0) != 0)
		{
		  mpz_add_ui(inv, inv, 1);
		  poly(inv, tmp);
		  mpz_tdiv_r_ui(tmp, tmp, factor_base[i].p);
		}
	      mpz_set(tmp, inv);
	    }
	  else
	    {
	      mpz_mul(tmp, tmp, inv);
	    }
	  mpz_tdiv_r_ui(tmp, tmp, factor_base[i].p);

	  if (j == 0)
	    factor_base[i].y = mpz_get_si(tmp);
	  else
	    factor_base[i].z = mpz_get_si(tmp);
	}

      /* Make sure that the two solutions are always positive */
      if (factor_base[i].y < 0)
	{
	  c = ((factor_base[i].y % factor_base[i].p) == 0 ? 0 : (-factor_base[i].y / factor_base[i].p) + 1);
	  factor_base[i].y = factor_base[i].y + c * factor_base[i].p;
	}
      if (factor_base[i].z < 0)
	{
	  c = ((factor_base[i].z % factor_base[i].p) == 0 ? 0 : (-factor_base[i].z / factor_base[i].p) + 1);
	  factor_base[i].z = factor_base[i].z + c * factor_base[i].p;
	}
    }

  mpz_clear(p); mpz_clear(tmp); mpz_clear(inv);
}

/* Compute the factor base */
void compute_factor_base()
{
  int i, n;
  double cp;
  mpz_t p, tmp;

  /* Initialize p to the first candidate for the factor base */
  /* We don't want two in the factor base since the Shanks-Tonelli algorithm can't handle an even prime */
  mpz_init(tmp);
  mpz_init_set_ui(p, 3);

  /* Find a factor base of base_size */
  i = 1;

  while (i < factor_base_size)
    {
      /* Only accept prime p if n is a quadratic residue modulo p (the Legendre symbol (n|p) = 1) */
      if (mpz_legendre(N, p) == 1)
        {
	  factor_base[i].p = mpz_get_ui(p);
          factor_base[i].logp = log(factor_base[i].p) / log(2);

	  /* Compute the square root of N mod p */
	  n = mpz_tdiv_r_ui(tmp, N, factor_base[i].p);
	  factor_base[i].sqrtn1 = tonelli(n, factor_base[i].p);
	  factor_base[i].sqrtn2 = factor_base[i].p - factor_base[i].sqrtn1;

	  i++;
        }
      mpz_nextprime(p, p);
    }

  /* Calculate the sieving threshold as log( M*Sqrt(N/2) / p_max^T ) */
  mpz_cdiv_q_ui(tmp, N, 2);
  mpz_sqrt(tmp, tmp);
  mpz_mul(tmp, tmp, M);
  cp = T * (log(factor_base[factor_base_size - 1].p) / log(2));
  mpz_ui_pow_ui(p, 2, (int)cp);
  mpz_cdiv_q(tmp, tmp, p);
  threshold = mpz_sizeinbase(tmp, 2);

  mpz_clear(tmp); mpz_clear(p);
}

/* Sieve */
int sieve()
{
  mpz_t f, q, r, start;
  int i, j, xindex_1, xindex_2, found, current_prime, current_prime_log, sieve_length, first, index;
  MPI_Status status;

  mpz_init(f); mpz_init(q); mpz_init(r); mpz_init(start);
  xindex_1 = xindex_2 = 0;

  /* Initialize the sieve subinterval */
  mpz_neg(start, M);
  mpz_sub(r, M, start);
  if (mpz_cmp_ui(r, SIEVE_LENGTH) > 0)
    sieve_length = SIEVE_LENGTH;
  else
    sieve_length = mpz_get_ui(r);

  /* Sieve the entire interval in steps on SIEVE_LENGTH */
  while (sieve_length > 0) {

    /* Initialize the sieve array with 0 */
    memset(logfxs, 0, sieve_length);

    /* Loop through all the elements in the factor base except the first one which is -1 */
    for (i = 1; i < factor_base_size; i++)
      {
	current_prime = factor_base[i].p;
	current_prime_log = factor_base[i].logp;

	for (j = 0; j < 2; j++)
	  {
	    /* Determine first x in the sieve interval that p divides */
	    first = (j == 0 ? factor_base[i].y : factor_base[i].z);
	    mpz_sub_ui(r, start, first);
	    mpz_tdiv_q_ui(q, r, current_prime);
	    mpz_mul_ui(q, q, current_prime);
	    mpz_add_ui(q, q, first);
	    mpz_sub(q, q, start);

	    if (j == 0)
	      xindex_1 = mpz_get_ui(q);
	    else
	      xindex_2 = mpz_get_ui(q);
	  }

	/* xindex_1 must be bigger than or equal to xindex_2 */
	if (xindex_1 < xindex_2)
	  {
	    int temp = xindex_2;
	    xindex_2 = xindex_1;
	    xindex_1 = temp;
	  }

	/* Loop through all x values in our interval - This is the real sieving */
	while (xindex_1 < sieve_length)
	  {
	    logfxs[xindex_1] = logfxs[xindex_1] + current_prime_log;
	    logfxs[xindex_2] = logfxs[xindex_2] + current_prime_log;
	    xindex_1 = xindex_1 + current_prime;
	    xindex_2 = xindex_2 + current_prime;
	  }
      }

    for (i = 0; i < sieve_length; i++)
      {
	/* If we are above the threshold do trial division with primes from the factor base */
	if (logfxs[i] >= threshold)
	  {
	    /* Calculate f(x) */
	    mpz_add_ui(r, start, i);
	    poly(r, f);

	    /* Clear the current row in the exponent matrix and take care of negative numbers */
	    mpz_set_ui(exponent_vector[0], 0);

	    if (mpz_cmp_ui(f, 0) < 0)
	      {
		mpz_neg(f, f);
		mpz_setbit(exponent_vector[0], 0);
	      }

	    /* Perform trial division and build an exponent vector */
	    found = (mpz_cmp_ui(f, 0) != 0);
	    while (found)
	      {
		found = 0;
		for (j = 1; j < factor_base_size; j++)
		  {
		    /* Check if the current element in the factor base divides f(x) */
		    if (mpz_cdiv_q_ui(q, f, factor_base[j].p) == 0)
		      {
			if (mpz_tstbit(exponent_vector[0], j))
			  mpz_clrbit(exponent_vector[0], j);
			else
			  mpz_setbit(exponent_vector[0], j);
			mpz_set(f, q);
			found = 1;
			break;
		      }
		  }
	      }

	    if (mpz_cmp_ui(f, 1) == 0)
	      {
		char buffer[MPI_BUFFER_SIZE];
		/* Wooohooooo, f(x) factors completely over the factor base */

	   	/* Send the result vector */
		send_mpz(RESULTTAG, 0, exponent_vector[0], buffer);

		/* Receive the master index */
		MPI_Recv(&index, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		/* Terminate if we got a DIETAG */
		if (status.MPI_TAG != RESULTACKTAG)
		  {
		    /* Clean up and return */
		    mpz_clear(f); mpz_clear(q); mpz_clear(r); mpz_clear(start);
		    return 1;
		  }

		/* Send ax + b to the master */
		mpz_add_ui(f, start, i);
		mpz_mul(f, f, A);
		mpz_add(f, f, B);
		send_mpz((index << 2) + AXBTAG, 0, f, buffer);

		/* Send a * f(x) to the master */
		mpz_add_ui(r, start, i);
		poly(r, f);
		mpz_mul(f, f, A);
		send_mpz((index << 2) + AFXTAG, 0, f, buffer);
	      }
	  }
      }

    /* Increase start and set sieve_length */
    mpz_add_ui(start, start, sieve_length);
    mpz_sub(r, M, start);
    if (mpz_cmp_ui(r, SIEVE_LENGTH) > 0)
      sieve_length = SIEVE_LENGTH;
    else
      sieve_length = mpz_get_ui(r);
  }

  /* Clean up */
  mpz_clear(f); mpz_clear(q); mpz_clear(r); mpz_clear(start);

  /* Return success */
  return 0;
}

/* Analyze the solutions in the matrix and print results */
void output_result()
{
  int first, i, j, retries, success;
  mpz_t x, y, p, q, f, gcd;

  mpz_init(x); mpz_init(y); mpz_init(p); mpz_init(q); mpz_init(gcd); mpz_init(f);

  /* Find first solution row in the matrix */
  mpz_set_ui(y, 2);
  mpz_pow_ui(y, y, factor_base_size);
  mpz_sub_ui(y, y, 1);
  first = retries = success = 0;
  for (i = factor_base_size + EXTRA - 1; i > 0; i--)
    {
      mpz_and(x, exponent_matrix[i], y);
      if (mpz_cmp_ui(x, 0) != 0)
	{
	  first = i + 1;
	  break;
	}
    }

  /* Loop through all possible solutions */
  for (i = first; i < factor_base_size + EXTRA; i++)
    {
      mpz_set_ui(y, 1);
      mpz_set_ui(x, 1);
      for (j = 0; j < factor_base_size + EXTRA; j++)
	{
	  if (mpz_tstbit(exponent_matrix[i], j + factor_base_size))
	    {
	      mpz_mul(y, y, afx[j]);
	      mpz_mul(x, x, axb[j]);
	    }
	}

      mpz_sqrt(y, y);
      mpz_sub(gcd, y, x);
      mpz_gcd(p, gcd, N);

      if (mpz_cmp_ui(p, 1) == 0 || mpz_cmp(p, N) == 0)
	/* We got the trivial solution. Let's try again */
	retries++;
      else {
	/* We found the non-trivial factorization */
	mpz_add(gcd, y, x);
	mpz_gcd(q, gcd, N);
	success = 1;
	break;
      }
    }

  /* Print the result */
  if (!success)
    {
      printf("Unable to find a factorization of ");
      mpz_out_str(stdout, 10, N);
      printf("\n\n");
    } else {
      if (retries > 0)
	printf("The first %d of %d solutions returned the trivial factorization\n", retries, factor_base_size + EXTRA - first);
      else
	printf("The first solution tried returned the non-trivial factorization\n");
      printf("The factorization is: \n\n");
      mpz_out_str(stdout, 10, N);
      printf(" = ");
      mpz_out_str(stdout, 10, p);
      printf(" * ");
      mpz_out_str(stdout, 10, q);
      printf("\n\n");
    }

  /* Clean up */
  mpz_clear(x); mpz_clear(y); mpz_clear(p); mpz_clear(q); mpz_clear(gcd); mpz_clear(f);
}

/* Handle abnormal termination of the master program */
static void catch_sigs(int signo)
{
  if (signo == SIGTERM)
    {
      printf("\nKilled...");
      printf("\e[?25h\n");
      MPI_Finalize();
      exit(0);
    }
  else if (signo == SIGINT)
    {
      printf("\e[?25h\n");
      MPI_Finalize();
      exit(0);
    }
}

void master(int argc, char **argv)
{
  time_t t_start, t_end, total_start, total_end;
  int relations, next_index, client_index, msg, length, nodes, terminated, done;
  char buffer[MPI_BUFFER_SIZE];
  MPI_Status status;
  struct sigaction sa_new;
  mpz_t tmp;

  /* Set up signal handling */
  sa_new.sa_handler = catch_sigs;
  sa_new.sa_flags = 0;
  sigemptyset(&sa_new.sa_mask);
  sigaddset(&sa_new.sa_mask, SIGTERM);
  sigaddset(&sa_new.sa_mask, SIGINT);
  sigaction(SIGINT, &sa_new, 0);
  sigaction(SIGTERM, &sa_new, 0);

  /* Turn off the cursor */
  printf("\e[?25l");

  /* Initialize variables */
  terminated = done = relations = next_index = 0;
  total_start = time(0);
  MPI_Comm_size(MPI_COMM_WORLD, &nodes);

  if (argc != 2)
    {
      printf("Usage: qs <number to factor>\n");
      done = 1;
    }
  else
    {
      printf("\n");
      length = strlen(argv[1]);
      if (!init_master(argv[1], length))
	{
	  printf("There are better ways to factor such a small number\n");
	  done = 1;
	}
      else
	{
	  /* Output banner greeting */
	  printf("Trying to factor ");
	  mpz_out_str(stdout, 10, N);
	  printf("\n\n");
	}
    }

  if (done == 0)
    {
      /* Check n for primality */
      if (mpz_probab_prime_p(N, 10) > 0)
	{
	  mpz_out_str(stdout, 10, N);
	  printf(" is a prime number\n\n");
	  done = 1;
	}
      else
	{
	  mpz_init(tmp);
	  t_start = time(0);

	  printf("Waiting for slave nodes to request work...\n\n");

	  printf("Sieving with %d slave nodes...\n", nodes-1);
	  printf("  Relations gathered: %5d, Needed: %5d\r", 0, factor_base_size + EXTRA);
	  fflush(stdout);
	  while (done == 0)
	    {
	      MPI_Recv(buffer, MPI_BUFFER_SIZE, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	      msg = status.MPI_TAG;

	      /* Receive ANY message */
	      if (msg == REQNTAG)
		{
		  if (next_index >= factor_base_size + EXTRA)
		    {
		      MPI_Send(buffer, 1, MPI_BYTE, status.MPI_SOURCE, DIETAG, MPI_COMM_WORLD);
		      terminated++;
		    } else {
		      /* Send N */
		      MPI_Send(&length, 1, MPI_INT, status.MPI_SOURCE, NLENGTHTAG, MPI_COMM_WORLD);
		      send_mpz(NTAG, status.MPI_SOURCE, N, buffer);
		    }
		}
	      else if (msg == REQDTAG)
		{
		  if (next_index >= factor_base_size + EXTRA)
                    {
		      MPI_Send(buffer, 1, MPI_BYTE, status.MPI_SOURCE, DIETAG, MPI_COMM_WORLD);
		      terminated++;
		    }
		  else
		    {
		      get_next_D();
		      /* Send D */
		      send_mpz(DTAG, status.MPI_SOURCE, D, buffer);
		    }
		}
	      else if (msg == RESULTTAG)
		{
		  /* Ask the client to terminate if we don't need more relations */
		  if (next_index >= factor_base_size + EXTRA)
		    {
		      MPI_Send(&client_index, 1, MPI_INT, status.MPI_SOURCE, DIETAG, MPI_COMM_WORLD);
		      terminated++;
		    }
		  else
		    {
		      /* Set next_index row in the matrix to the value received */
		      buffer2mpz(tmp, buffer);
		      mpz_set(exponent_matrix[next_index], tmp);

		      /* Send a reply to the client containing the index of the next two results */
		      client_index = next_index;
		      MPI_Send(&client_index, 1, MPI_INT, status.MPI_SOURCE, RESULTACKTAG, MPI_COMM_WORLD);
		      next_index++;
		    }
		}

	      client_index = msg >> 2;
	      msg = msg & 3;

	      switch (msg)
		{
		case AXBTAG:
		  /* Set client_index entry in the axb array to the value received */
		  buffer2mpz(axb[client_index], buffer);

		  break;
		case AFXTAG:
		  /* Set client_index entry in the afx array to the value received */
		  buffer2mpz(afx[client_index], buffer);

		  relations++;
		  done = (relations >= factor_base_size + EXTRA);
		  printf("  Relations gathered: %5d, Needed: %5d\r", relations, factor_base_size + EXTRA);
		  fflush(stdout);
		  break;
		}
	    }

	  mpz_clear(tmp);
	  t_end = time(0);
	  printf("\nDone (%d seconds)\n\n", (int)(t_end - t_start));

	  /* Gauss elimination */
	  printf("Performing gauss elimination...\n");
	  t_start = time(0);
	  gauss_eliminate(factor_base_size + EXTRA, exponent_matrix, EXTRA);
	  t_end = time(0);
	  printf("Done (%d seconds)\n\n", (int)(t_end - t_start));

	  /* Output result */
	  output_result();
	}

      total_end = time(0);
      printf("Total running time: %d seconds\n", (int)(total_end - total_start));
    }

  printf("\nWaiting for clients to terminate...\n");
  printf("  %d of %d clients have terminated\r", terminated, nodes-1);
  fflush(stdout);
  while (terminated != nodes - 1)
    {
      MPI_Recv(buffer, MPI_BUFFER_SIZE, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      MPI_Send(buffer, 1, MPI_BYTE, status.MPI_SOURCE, DIETAG, MPI_COMM_WORLD);
      terminated++;
      printf("  %d of %d clients have terminated\r", terminated, nodes-1);
      fflush(stdout);
    }
  printf("\n");

  return;
}

void slave()
{
  int msg, length;
  char buffer[MPI_BUFFER_SIZE];
  MPI_Status status;

  /* Get n */
  MPI_Send(&msg, 1, MPI_INT, 0, REQNTAG, MPI_COMM_WORLD);

  MPI_Recv(&length, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  if (status.MPI_TAG != NLENGTHTAG) return;
  MPI_Recv(buffer, MPI_BUFFER_SIZE, MPI_BYTE, 0, NTAG, MPI_COMM_WORLD, &status);

  /* Get n from the message */
  buffer2mpz(N, buffer);

  init_slave(length);
  compute_factor_base();

  /* Get D */
  MPI_Send(&msg, 1, MPI_INT, 0, REQDTAG, MPI_COMM_WORLD);
  MPI_Recv(buffer, MPI_BUFFER_SIZE, MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
  if (status.MPI_TAG != DTAG) return;

  /* Get D from the message */
  buffer2mpz(D, buffer);

  while(1)
    {
      generate_new_polynomial();
      initialize_polynomial();
      if (sieve() != 0) return;

      /* Request new D */
      MPI_Send(&msg, 1, MPI_INT, 0, REQDTAG, MPI_COMM_WORLD);

      MPI_Recv(buffer, MPI_BUFFER_SIZE, MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      if (status.MPI_TAG != DTAG) return;

      /* Get D from the message */
      buffer2mpz(D, buffer);
    }
}

int main(int argc, char **argv)
{
  /*

    Test Numbers taken from Olof Aasbrink and Joel Brynielsson's paper about the QS

    T20: 18567078082619935259
    T30: 350243405507562291174415825999
    T40: 5705979550618670446308578858542675373983
    T45: 732197471686198597184965476425281169401188191
    T50: 53468946676763197941455249471721044636943883361749
    T55: 5945326581537513157038636316967257854322393895035230547
    T60: 676292275716558246502605230897191366469551764092181362779759

    This number is from the RSA challenge and currently unfactored. Knock yourself out

    RSA-576: 188198812920607963838697239461650439807163563379417
             382700763356422988859715234665485319060606504743045
	     317388011303396716199692321205734031879550656996221
	     3051687593076502570

  */

  int myrank;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == 0)
    {
      master(argc, argv);
      catch_sigs(SIGINT);
    }
  else
    {
      slave();
    }

  MPI_Finalize();
  return 0;
}
