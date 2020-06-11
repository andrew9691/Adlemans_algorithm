#ifndef GAUSS_H
#define GAUSS_H

#include <ostream>
#include <vector>
#include <NTL\ZZ.h>

std::vector<NTL::ZZ> myGauss(std::vector<std::vector<NTL::ZZ>> mat, std::vector<NTL::ZZ> r, const NTL::ZZ& p)
{
	std::vector<NTL::ZZ> res;
	size_t k = 0;
	size_t n = r.size();
	NTL::ZZ max;
	int index;
	while (k < n)
	{
		max = NTL::abs(mat[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (NTL::abs(mat[i][k]) > max)
			{
				max = NTL::abs(mat[i][k]);
				index = i;
			}
		}

		if (max == (NTL::ZZ)0)
			return res;

		for (int j = 0; j < n; j++)
		{
			NTL::ZZ temp = mat[k][j];
			mat[k][j] = mat[index][j];
			mat[index][j] = temp;
		}
		NTL::ZZ temp = r[k];
		r[k] = r[index];
		r[index] = temp;

		for (int i = 0; i < n; i++)
		{
			NTL::ZZ temp = mat[i][k];
			if (temp == (NTL::ZZ)0)
				continue;

			if (i == k)
				continue;

			for (int j = 0; j < n; j++)
				mat[i][j] = mat[k][k] * mat[i][j] - temp * mat[k][j];
			r[i] = mat[k][k] * r[i] - temp * r[k];
		}
		k++;
	}

	for (int i = 0; i < n; i++)
	{
		NTL::ZZ gcd = NTL::GCD(mat[i][i], r[i]);
		mat[i][i] /= gcd;
		r[i] /= gcd;

		if (mat[i][i] < 0)
		{
			r[i] *= -1;
			mat[i][i] *= -1;
		}
		res.emplace_back(r[i] * NTL::InvMod(mat[i][i], p - 1));
	}

	return res;
}

#endif
