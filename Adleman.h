#ifndef ADLEMAN_H
#define ADLEMAN_H

#include <NTL\ZZ_p.h>
#include <pair.h>
#include <set>
#include <xutility>
#include "Gauss.h"

std::vector<std::pair<NTL::ZZ, NTL::ZZ>> B_smooth_factorization(const NTL::ZZ& x, const std::vector<NTL::ZZ>& factor_base)
{
	NTL::ZZ a = x;
	NTL::ZZ prime_div = (NTL::ZZ)2;
	NTL::ZZ deg = (NTL::ZZ)0;

	std::vector<std::pair<NTL::ZZ, NTL::ZZ>> res;
	while (a > 1)
	{
		if (prime_div > factor_base[factor_base.size() - 1])
		{
			res.erase(res.begin(), res.end());
			return res;
		}

		while (a % prime_div == 0)
		{
			a /= prime_div;
			deg++;
		}

		res.emplace_back(std::make_pair(prime_div, deg));

		deg = 0;
		prime_div = NTL::NextPrime(prime_div + 1);
	}

	if (res.size() < factor_base.size())
	{
		while (res[res.size() - 1].first != factor_base[factor_base.size() - 1])
			res.emplace_back(std::make_pair(NTL::NextPrime(res[res.size() - 1].first + 1), (NTL::ZZ)0));
	}

	return res;
}

bool equal_congruences(const std::vector<std::pair<NTL::ZZ, NTL::ZZ>>& c1, const std::vector<std::pair<NTL::ZZ, NTL::ZZ>>& c2)
{
	std::vector<NTL::ZZ> v1, v2;
	for (int i = 0; i < c1.size(); i++)
	{
		if (c1[i].second > 0)
			v1.emplace_back(c1[i].first);

		if (c2[i].second > 0)
			v2.emplace_back(c2[i].first);
	}

	if (v1.size() != v2.size())
		return false;

	for (int i = 0; i < v1.size(); i++)
	{
		if (v1[i] != v2[i])
			return false;
	}
	return true;
}

NTL::ZZ adleman_alg(const NTL::ZZ& a, const NTL::ZZ& b, const NTL::ZZ& p)
{
	NTL::ZZ_p::init(p);

	NTL::ZZ B = NTL::conv<NTL::ZZ>(std::round(std::pow(2, std::sqrt(std::log2(NTL::conv<double>(p))))));

	std::set<NTL::ZZ> used_factors;
	std::vector<NTL::ZZ> factor_base;
	NTL::ZZ prime = (NTL::ZZ)2;
	while (prime <= B)
	{
		factor_base.emplace_back(prime);
		prime = NTL::NextPrime(prime + 1);
	}

	std::vector<std::vector<std::pair<NTL::ZZ, NTL::ZZ>>> factored_apowr;
	std::vector<NTL::ZZ_p> r_i;

	NTL::ZZ_p r = (NTL::ZZ_p)1;
	int flag = 0;
	while (r_i.size() != factor_base.size())
	{
		m:
		r *= NTL::random_ZZ_p();
		if (r == NTL::conv<NTL::ZZ_p>(p - 1) || r == 0)
			goto m;

		std::vector<std::pair<NTL::ZZ, NTL::ZZ>> factored_r = B_smooth_factorization(NTL::PowerMod(a, NTL::conv<NTL::ZZ>(r), p), factor_base);
		if (factored_r.size() > 0)
		{
			if (flag == 0)
			{
				int f2 = 0;
				for (int i = 0; i < factored_r.size(); i++)
				{
					if (factored_r[i].second > 0 && std::find(used_factors.begin(), used_factors.end(), factored_r[i].first) == used_factors.end())
					{
						f2 = 1;
						used_factors.insert(factored_r[i].first);
					}
				}

				if (used_factors.size() == factor_base.size())
					flag = 1;

				if (f2 == 1)
				{
					r_i.emplace_back(r);
					factored_apowr.emplace_back(factored_r);
				}
			}
			else
			{
				int f1 = 1;
				for (int i = 0; i < factored_apowr.size(); i++)
				{
					if (equal_congruences(factored_apowr[i], factored_r))
					{
						f1 = 0;
						break;
					}
				}

				if (f1 == 1)
				{
					r_i.emplace_back(r);
					factored_apowr.emplace_back(factored_r);
				}
			}
		}
	}

	int n = r_i.size();
	std::vector<std::vector<NTL::ZZ>> mat;
	std::vector<NTL::ZZ> sub_mat;
	std::vector<NTL::ZZ> y;
	for (int i = 0; i < n; i++)
	{
		y.emplace_back(NTL::conv<NTL::ZZ>(r_i[i]));

		for (int j = 0; j < n; j++)
			sub_mat.emplace_back(factored_apowr[i][j].second);

		mat.emplace_back(sub_mat);
		sub_mat.erase(sub_mat.begin(), sub_mat.end());
	}
	std::vector<NTL::ZZ> X = myGauss(mat, y, p);

	NTL::ZZ prmdv = (NTL::ZZ)2;
	for (int i = 0; i < n; i++)
		prmdv = NTL::NextPrime(prmdv + 1);

	NTL::ZZ q;
	std::vector<std::pair<NTL::ZZ, NTL::ZZ>> factored_apowr_b;

	while (factored_apowr_b.size() == 0)
	{
		m4:
		NTL::ZZ_p r = NTL::random_ZZ_p();
		if (r == NTL::conv<NTL::ZZ_p>(p - 1) || r == 0)
			goto m4;

		NTL::ZZ rii = NTL::conv<NTL::ZZ>(r);
		NTL::ZZ powa = NTL::PowerMod(a, rii, p);
		NTL::ZZ_p apowr_b = NTL::conv<NTL::ZZ_p>(powa * b);

		q = NTL::conv<NTL::ZZ>(r);
		factored_apowr_b = B_smooth_factorization(NTL::conv<NTL::ZZ>(apowr_b), factor_base);
	}

	NTL::ZZ_p::init(p - 1);
	NTL::ZZ_p res = NTL::conv<NTL::ZZ_p>(-q);

	for (int i = 0; i < factored_apowr_b.size(); i++)
		res += NTL::conv<NTL::ZZ_p>(factored_apowr_b[i].second * X[i]);

	return NTL::conv<NTL::ZZ>(res);
}

#endif
