#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>

// represents n/d
typedef struct {
	int64_t n;
	int64_t d;	// >0
} rational_t;

int64_t gcd_r(int64_t a, int64_t b)
{
	return (b == 0) ? a : gcd_r(b, a % b);
}

int64_t gcd(int64_t a, int64_t b)
{
	return (a > b) ? gcd_r(a, b) : gcd_r(b, a);
}

int64_t lcm(int64_t a, int64_t b)
{
	return (a*b)/gcd(a, b);
}

int64_t abs64(int64_t a)
{
	return (a < 0) ? (-a) : a;
}

void rational_assign(rational_t *r, int64_t n, int64_t d)
{
	if (d < 0) {
		n = -n;
		d = -d;
	}

	int64_t const g = gcd(abs64(n), d);
	r->n = n/g;
	r->d = d/g;
}

rational_t make_rational(int64_t n, int64_t d)
{
	rational_t r;
	rational_assign(&r, n, d);
	return r;
}

void rational_mul(rational_t *r, rational_t a)
{
	rational_assign(r, r->n*a.n, r->d*a.d);
}

void rational_add(rational_t *r, rational_t a)
{
	rational_assign(r, r->n*a.d + a.n*r->d, r->d*a.d);
}

int64_t largest_square_factor(int64_t *n)
{
	int64_t const n_orig = *n;
	int64_t remain = n_orig;
	int64_t factor = 1;
	for (int64_t i = 2; i*i <= remain; ++i) {
		int64_t const ii = i*i;
		while ((remain % ii) == 0) {
			factor *= i;
			remain /= ii;
		}
	}
	if (n_orig < 0) {
		remain = -remain;
	}
	*n = remain;
	return factor;
}

void rational_extract_square(rational_t *a2, rational_t *b)
{
	rational_t const b_orig = *b;
	rational_mul(a2, make_rational(b_orig.n*b_orig.n, b_orig.d*b_orig.d));
	int64_t n = largest_square_factor(&a2->n);
	if (b_orig.n < 0) {
		n = -n;
	}
	int64_t const d = largest_square_factor(&a2->d);
	*b = make_rational(n, d);
}

void rational_tex(rational_t r, bool is_first, bool is_constant)
{
	if (r.n < 0) {
		printf("-");
	} else if (!is_first) {
		printf("+");
	}

	int64_t const abs_n = abs64(r.n);
	if (r.d == 1) {
		if (abs_n == 1 && !is_constant) {
			// omit
		} else {
			printf("%"PRId64, abs_n);
		}
	} else {
		printf("\\frac{%"PRId64"}{%"PRId64"}", abs_n, r.d);
	}
}

// represents some power of 2 variables
typedef struct {
	int p0;
	int p1;
} var_t;

var_t make_var(int p0, int p1)
{
	return (var_t){ .p0 = p0, .p1 = p1 };
}

void var_mul(var_t *v, var_t a)
{
	v->p0 += a.p0;
	v->p1 += a.p1;
}

bool var_is_equal(var_t a, var_t b)
{
	return (a.p0 == b.p0) && (a.p1 == b.p1);
}

bool var_is_constant(var_t v)
{
	return (v.p0 == 0 && v.p1 == 0);
}

void var_tex(
	var_t v,
	char const *pre0,
	char const *pre1,
	char const *post)
{
	if (v.p0 == 1) {
		printf("%s%s", pre0, post);
	} else if (v.p0 != 0) {
		printf("%s^%d%s", pre0, v.p0, post);
	}

	if (v.p0 != 0 && v.p1 != 0) {
		printf("\\,");
	}

	if (v.p1 == 1) {
		printf("%s%s", pre1, post);
	} else if (v.p1 != 0) {
		printf("%s^%d%s", pre1, v.p1, post);
	}
}

// some coefficient of a polynomial in two variables
typedef struct
{
	rational_t r;
	var_t v;
} coeff_t;

// polynomial in two variables
typedef struct
{
	coeff_t *coeffs;
	int size;
} polynomial_t;

void polynomial_create_rational(polynomial_t *p, rational_t r)
{
	p->coeffs = (coeff_t *)malloc(sizeof(coeff_t));
	p->coeffs[0] = (coeff_t){ .r = r, .v = { 0, 0 } };
	p->size = 1;
}

void polynomial_create_coeff(polynomial_t *p, rational_t r, var_t v)
{
	p->coeffs = (coeff_t *)malloc(sizeof(coeff_t));
	p->coeffs[0] = (coeff_t){ .r = r, .v = v };
	p->size = 1;
}

void polynomial_destroy(polynomial_t *p)
{
	free(p->coeffs);
	p->coeffs = NULL;
	p->size = 0;
}

void polynomial_mul_rational(polynomial_t *p, rational_t r)
{
	for (int i = 0, sz = p->size; i < sz; ++i) {
		rational_mul(&p->coeffs[i].r, r);
	}
}

void polynomial_mul_var(polynomial_t *p, var_t v)
{
	for (int i = 0, sz = p->size; i < sz; ++i) {
		var_mul(&p->coeffs[i].v, v);
	}
}

void polynomial_mul_coeff(polynomial_t *p, rational_t r, var_t v)
{
	for (int i = 0, sz = p->size; i < sz; ++i) {
		rational_mul(&p->coeffs[i].r, r);
		var_mul(&p->coeffs[i].v, v);
	}
}

void polynomial_add_coeff(polynomial_t *p, rational_t r, var_t v)
{
	for (int i = 0, sz = p->size; i < sz; ++i) {
		if (var_is_equal(p->coeffs[i].v, v)) {
			rational_add(&p->coeffs[i].r, r);
			return;
		}
	}
	int const old_sz = p->size;
	int const new_sz = old_sz + 1;
	p->coeffs = realloc(p->coeffs, new_sz*sizeof(coeff_t));
	p->coeffs[old_sz] = (coeff_t){ .r = r, .v = v };
	p->size = new_sz;
}

void polynomial_add_polynomial(polynomial_t *p, polynomial_t const *q)
{
	for (int i = 0, sz = q->size; i < sz; ++i) {
		polynomial_add_coeff(p, q->coeffs[i].r, q->coeffs[i].v);
	}
}

void polynomial_extract_factor(polynomial_t *p, rational_t *r)
{
	if (p->size == 1) {
		rational_mul(r, p->coeffs[0].r);
		if (var_is_constant(p->coeffs[0].v)) {
			p->size = 0;
		} else {
			p->coeffs[0].r = make_rational(1, 1);
		}
	} else {
		int64_t n_gcd = 0;
		int64_t d_lcm = 1;
		for (int i = 0, sz = p->size; i < sz; ++i) {
			rational_t const t = p->coeffs[i].r;
			n_gcd = gcd(n_gcd, abs64(t.n));
			d_lcm = lcm(d_lcm, t.d);
		}
		rational_mul(r, make_rational(n_gcd, d_lcm));
		for (int i = 0, sz = p->size; i < sz; ++i) {
			rational_mul(&p->coeffs[i].r, make_rational(d_lcm, n_gcd));
		}
	}
}

void polynomial_tex(
	polynomial_t const *p,
	char const *pre0,
	char const *pre1,
	char const *post)
{
	if (p->size > 1) {
		printf("(");
	}
	for (int i = 0, sz = p->size; i < sz; ++i) {
		rational_tex(p->coeffs[i].r, i == 0, var_is_constant(p->coeffs[i].v));
		var_tex(p->coeffs[i].v, pre0, pre1, post);
	}
	if (p->size > 1) {
		printf(")");
	}
}

#define make_var_sin() make_var(1, 0)
#define make_var_cos() make_var(0, 1)

void legendre_p(polynomial_t *p, int l, int m)
{
	if (l == 0 && m == 0) {
		polynomial_create_rational(p, make_rational(1, 1));
	} else if (l == m) {
		legendre_p(p, m - 1, m - 1);
		polynomial_mul_coeff(p, make_rational(2*m - 1, 1), make_var_sin());
	} else if (l == m + 1) {
		legendre_p(p, m, m);
		polynomial_mul_coeff(p, make_rational(2*m + 1, 1), make_var_cos());
	} else {
		legendre_p(p, l - 1, m);
		polynomial_mul_coeff(p, make_rational(2*l - 1, l - m), make_var_cos());

		polynomial_t q;
		legendre_p(&q, l - 2, m);
		polynomial_mul_rational(&q, make_rational(1 - l - m, l - m));

		polynomial_add_polynomial(p, &q);
		polynomial_destroy(&q);
	}
}

void sin_p(polynomial_t *p, int m);
void cos_p(polynomial_t *p, int m);

void phi_p(polynomial_t *p, int m)
{
	if (m == 0) {
		polynomial_create_rational(p, make_rational(1, 1));
	} else if (m < 0) {
		sin_p(p, -m);
	} else {
		cos_p(p, m);
	}
}

void sin_p(polynomial_t *p, int m)
{
	if (m == 1) {
		polynomial_create_coeff(p, make_rational(1, 1), make_var_sin());
	} else {
		cos_p(p, m - 1);
		polynomial_mul_var(p, make_var_sin());

		polynomial_t q;
		sin_p(&q, m - 1);
		polynomial_mul_var(&q, make_var_cos());

		polynomial_add_polynomial(p, &q);
		polynomial_destroy(&q);
	}
}

void cos_p(polynomial_t *p, int m)
{
	if (m == 1) {
		polynomial_create_coeff(p, make_rational(1, 1), make_var_cos());
	} else {
		cos_p(p, m - 1);
		polynomial_mul_var(p, make_var_cos());

		polynomial_t q;
		sin_p(&q, m - 1);
		polynomial_mul_coeff(&q, make_rational(-1, 1), make_var_sin());

		polynomial_add_polynomial(p, &q);
		polynomial_destroy(&q);
	}
}

int factorial(int a)
{
	return (a < 2) ? 1 : (a*factorial(a - 1));
}

int main(int argc, char *argv[])
{
	int max_order = 2;
	if (argc == 2) {
		max_order = atoi(argv[1]);
	}

	for (int l = 0; l <= max_order; ++l)
	for (int m = -l; m <= l; ++m) {
		int const abs_m = abs(m);

		polynomial_t p;
		legendre_p(&p, l, abs_m);

		polynomial_t q;
		phi_p(&q, m);

		rational_t n2 = make_rational(
			(2*l + 1)*factorial(l - abs_m)*((m == 0) ? 1 : 2),
			4*factorial(l + abs_m));

		rational_t k = make_rational(1, 1);
		polynomial_extract_factor(&p, &k);
		polynomial_extract_factor(&q, &k);
		rational_extract_square(&n2, &k);

		printf("Y^%d_{%d} &= ", l, m);

		rational_tex(k, true, false);

		printf("\\sqrt{\\frac{%"PRId64"}{", n2.n);
		if (n2.d != 1) {
			printf("%"PRId64, n2.d);
		}
		printf("\\pi}}");

		polynomial_tex(&q, "\\sin", "\\cos", "\\phi");
		if (q.size == 1 && p.size == 1) {
			printf("\\,");
		}
		polynomial_tex(&p, "\\sin", "\\cos", "\\theta");

		printf("\\\\\n");

		polynomial_destroy(&q);
		polynomial_destroy(&p);
	}
	return 0;
}
