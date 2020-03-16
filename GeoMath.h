
class MathsFun {	
public:

	float max(float x, float y) {
		return (x > y) ? x : y;
	}

	float min(float x, float y) {
		return (x < y) ? x : y;
	}

	template <typename T>
	T myMax(T x, T y) {
		return (x > y) ? x : y;
	}

	template <typename T>
	T myMin(T x, T y) {
		return (x < y) ? x : y;
	}

	
	template <typename T>
		void sortem(int ib, int ie, fvec& a, int iperm, vect<T>& b, vect<T>&c) {

		/* you can pass thesame in twice 
		-----------------------------------------------------------------------

							  Quickersort Subroutine
							  **********************

		 This is a subroutine for sorting a real array in ascending order. This
		 is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
		 in collected algorithms of the ACM.

		 The method used is that of continually splitting the array into parts
		 such that all elements of one part are less than all elements of the
		 other, with a third part in the middle consisting of one element.  An
		 element with value t is chosen arbitrarily (here we choose the middle
		 element). i and j give the lower and upper limits of the segment being
		 split.  After the split a value q will have been found such that
		 a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
		 performs operations on the two segments (i,q-1) and (q+1,j) as follows
		 The smaller segment is split and the position of the larger segment is
		 stored in the lt and ut arrays.  If the segment to be split contains
		 two or fewer elements, it is sorted and another segment is obtained
		 from the lt and ut arrays.  When no more segments remain, the array
		 is completely sorted.


		 INPUT PARAMETERS:

		   ib,ie        start and end index of the array to be sorteda
		   a            array, a portion of which has to be sorted.
		   iperm        0 no other array is permuted.
						1 array b is permuted according to array a
						2 arrays b,c are permuted.
						3 arrays b,c,d are permuted.
						4 arrays b,c,d,e are permuted.
						5 arrays b,c,d,e,f are permuted.
						6 arrays b,c,d,e,f,g are permuted.
						7 arrays b,c,d,e,f,g,h are permuted.
					   >7 no other array is permuted.

		   b,c,d,e,f,g,h  arrays to be permuted according to array a.

		 OUTPUT PARAMETERS:

			a              = the array, a portion of which has been sorted.

			b,c,d,e,f,g,h  = arrays permuted according to array a (see iperm)

		 NO EXTERNAL ROUTINES REQUIRED:

		%-----------------------------------------------------------------------
		*/

		fvec  d, e, f, g, h;
		d.zeros(1); f.zeros(1); g.zeros(1), h.zeros(1);;

		// The dimensions for lt and ut have to be at least log (base 2) n
		int iring, lt[64], ut[64], i, j, k, m, p, q;
		float ta, tb, tc, td, te, tf, tg, th, xa, xb, xc, xd, xe, xf, xg, xh;

		// Initialize:
		j = ie;
		m = 1;
		i = ib;
		iring = iperm + 1;
		if (iperm > 7) iring = 1;

		// If this segment has more than two elements  we split it
	A: if (j - i - 1 < 0)
		goto B;
	   else if (j - i - 1 == 0)
		goto D;
	   else if (j - i - 1 > 0)
		goto K;

	   // p is the position of an arbitrary element in the segment we choose the
	   // middle element. Under certain circumstances it may be advantageous
	   // to choose p at random.
   K: p = (j + i) / 2;
	   ta = a[p];
	   a[p] = a[i];
	   if (iring >= 2) {
		   tb = b[p];
		   b[p] = b[i];
	   }
	   if (iring >= 3) {
		   tc = c[p];
		   c[p] = c[i];
	   }
	   if (iring >= 4) {
		   td = d[p];
		   d[p] = d[i];
	   }
	   if (iring >= 5) {
		   te = e[p];
		   e[p] = e[i];
	   }
	   if (iring >= 6) {
		   tf = f[p];
		   f[p] = f[i];
	   }
	   if (iring >= 7) {
		   tg = g[p];
		   g[p] = g[i];
	   }
	   if (iring >= 8) {
		   th = h[p];
		   h[p] = h[i];
	   }

	   //Start at the beginning of the segment, search for k such that a[k]>t
	   q = j;
	   k = i;
   H: k = k + 1;
	   if (k > q) goto F;
	   if (a[k] <= ta) goto H;

	   // Such an element has now been found now search for a q such that a[q]<t
	   // starting at the end of the segment.
   I:
	   if (a[q] < ta) goto J;
	   q = q - 1;
	   if (q > k)     goto I;
	   goto E;

	   // a[q] has now been found. we interchange a[q] and a[k]
   J: xa = a[k];
	   a[k] = a[q];
	   a[q] = xa;
	   if (iring >= 2) {
		   xb = b[k];
		   b[k] = b[q];
		   b[q] = xb;
	   }
	   if (iring >= 3) {
		   xc = c[k];
		   c[k] = c[q];
		   c[q] = xc;
	   }
	   if (iring >= 4) {
		   xd = d[k];
		   d[k] = d[q];
		   d[q] = xd;
	   }
	   if (iring >= 5) {
		   xe = e[k];
		   e[k] = e[q];
		   e[q] = xe;
	   }
	   if (iring >= 6) {
		   xf = f[k];
		   f[k] = f[q];
		   f[q] = xf;
	   }
	   if (iring >= 7) {
		   xg = g[k];
		   g[k] = g[q];
		   g[q] = xg;
	   }
	   if (iring >= 8) {
		   xh = h[k];
		   h[k] = h[q];
		   h[q] = xh;
	   }

	   // Update q and search for another pair to interchange:
	   q = q - 1;
	   goto H;
   E: q = k - 1;
   F:

	   // The upwards search has now met the downwards search:
	   a[i] = a[q];
	   a[q] = ta;
	   if (iring >= 2) {
		   b[i] = b[q];
		   b[q] = tb;
	   }
	   if (iring >= 3) {
		   c[i] = c[q];
		   c[q] = tc;
	   }
	   if (iring >= 4) {
		   d[i] = d[q];
		   d[q] = td;
	   }
	   if (iring >= 5) {
		   e[i] = e[q];
		   e[q] = te;
	   }
	   if (iring >= 6) {
		   f[i] = f[q];
		   f[q] = tf;
	   }
	   if (iring >= 7) {
		   g[i] = g[q];
		   g[q] = tg;
	   }
	   if (iring >= 8) {
		   h[i] = h[q];
		   h[q] = th;
	   }

	   // The segment is now divided in three parts: (i,q-1),[q],(q+1,j)
	   // store the position of the largest segment in lt and ut
	   if (2 * q <= i + j) goto G;
	   lt[m] = i;
	   ut[m] = q - 1;
	   i = q + 1;
	   goto C;
   G: lt[m] = q + 1;
	   ut[m] = j;
	   j = q - 1;

	   //Update m and split the new smaller segment
   C: m = m + 1;
	   goto A;

	   // We arrive here if the segment has  two elements we test to see if
	   // the segment is properly ordered if not, we perform an interchange
   D:
	   if (a[i] <= a[j]) goto B;
	   xa = a[i];
	   a[i] = a[j];
	   a[j] = xa;
	   if (iring >= 2) {
		   xb = b[i];
		   b[i] = b[j];
		   b[j] = xb;
	   }
	   if (iring >= 3) {
		   xc = c[i];
		   c[i] = c[j];
		   c[j] = xc;
	   }
	   if (iring >= 4) {
		   xd = d[i];
		   d[i] = d[j];
		   d[j] = xd;
	   }
	   if (iring >= 5) {
		   xe = e[i];
		   e[i] = e[j];
		   e[j] = xe;
	   }
	   if (iring >= 6) {
		   xf = f[i];
		   f[i] = f[j];
		   f[j] = xf;
	   }
	   if (iring >= 7) {
		   xg = g[i];
		   g[i] = g[j];
		   g[j] = xg;
	   }
	   if (iring >= 8) {
		   xh = h[i];
		   h[i] = h[j];
		   h[j] = xh;
	   }

	   // If lt and ut contain more segments to be sorted repeat process:
   B: m = m - 1;
	   if (m <= 0) return;
	   i = lt[m];
	   j = ut[m];
	   goto A;

	}



};



