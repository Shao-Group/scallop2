/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "util.h"

vector<int> get_random_permutation(int n)
{
	vector<int> v;
	for(int i = 0; i < n; i++) v.push_back(i);
	for(int i = 0; i < n; i++)
	{
		int k = rand() % (n - i);
		int x = v[k];
		v[k] = v[n - i - 1];
		v[n - i - 1] = x;
	}
	return v;
}

size_t string_hash(const std::string& str)
{
	size_t hash = 1315423911;
	for(std::size_t i = 0; i < str.length(); i++)
	{
		hash ^= ((hash << 5) + str[i] + (hash >> 2));
	}

	return (hash & 0x7FFFFFFF);
}

size_t vector_hash(const vector<int32_t> & vec)
{
	size_t seed = vec.size();
	for(int i = 0; i < vec.size(); i++)
	{
		seed ^= (size_t)(vec[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	}
	return (seed & 0x7FFFFFFF);
}
