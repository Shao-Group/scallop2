/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include <cassert>
#include "transcript_set.h"
#include "config.h"

trans_item::trans_item()
{}

trans_item::trans_item(const transcript &t, int c, int s)
{
	trst = t;
	count = c;
	if(samples.find(s) == samples.end()) samples.insert(make_pair(s, t.coverage));
	else if(samples[s] < t.coverage) samples[s] = t.coverage;
}

int trans_item::merge(const trans_item &ti, int mode)
{
	if(mode == TRANSCRIPT_COUNT_ADD_COVERAGE_ADD) 
	{
		trst.coverage += ti.trst.coverage;

		trst.extend_bounds(ti.trst);
		count += ti.count;

		for(auto &x : ti.samples)
		{
			if(samples.find(x.first) == samples.end()) samples.insert(x);
			else if(samples[x.first] < x.second) samples[x.first] = x.second;
		}
		//samples.insert(ti.samples.begin(), ti.samples.end());
	}
	else if(mode == TRANSCRIPT_COUNT_ADD_COVERAGE_MAX) 
	{
		if(trst.coverage < ti.trst.coverage) trst.coverage = ti.trst.coverage;

		trst.extend_bounds(ti.trst);
		count += ti.count;

		for(auto &x : ti.samples)
		{
			if(samples.find(x.first) == samples.end()) samples.insert(x);
			else if(samples[x.first] < x.second) samples[x.first] = x.second;
		}
		//samples.insert(ti.samples.begin(), ti.samples.end());
	}
	else if(mode == TRANSCRIPT_COUNT_ADD_COVERAGE_MIN) 
	{
		if(trst.coverage > ti.trst.coverage) trst.coverage = ti.trst.coverage;

		trst.extend_bounds(ti.trst);
		count += ti.count;

		for(auto &x : ti.samples)
		{
			if(samples.find(x.first) == samples.end()) samples.insert(x);
			else if(samples[x.first] < x.second) samples[x.first] = x.second;
		}
		//samples.insert(ti.samples.begin(), ti.samples.end());
	}
	else if(mode == TRANSCRIPT_COUNT_ADD_COVERAGE_NUL) 
	{
		count += ti.count;
	}
	else assert(false);

	return 0;
}

int merge_sorted_trans_items(vector<trans_item> &vx, const vector<trans_item> &vy, int mode1, int mode2, double single_exon_ratio)
{
	vector<trans_item> vz;
	int kx = 0, ky = 0;
	while(kx < vx.size() && ky < vy.size())
	{
		int b = vx[kx].trst.compare1(vy[ky].trst, single_exon_ratio);
		if(b == 0)
		{
			if(vx[kx].trst.exons.size() <= 1) vx[kx].merge(vy[ky], mode1);
			if(vx[kx].trst.exons.size() >= 2) vx[kx].merge(vy[ky], mode2);
			vz.push_back(vx[kx]);
			kx++;
			ky++;
		}
		else if(b == 1)
		{
			vz.push_back(vx[kx]);
			kx++;
		}
		else if(b == -1)
		{
			vz.push_back(vy[ky]);
			ky++;
		}
		else assert(false);
	}

	assert(kx == vx.size() || ky == vy.size());

	for(int i = kx; i < vx.size(); i++) vz.push_back(vx[i]);
	for(int i = ky; i < vy.size(); i++) vz.push_back(vy[i]);

	vx.clear();
	vx = vz;
	return 0;
}

transcript_set::transcript_set(const string &c, double s)
{
	chrm = c;
	single_exon_overlap = s;
}

transcript_set::transcript_set(const transcript &t, int count, int sid, double overlap)
{
	chrm = t.seqname;
	single_exon_overlap = overlap;

	size_t h = t.get_intron_chain_hashing();
	trans_item ti(t, count, sid);
	vector<trans_item> v;
	v.push_back(std::move(ti));

	mt.insert(make_pair(h, v));
}

int transcript_set::add(const transcript &t, int count, int sid, int mode1, int mode2)
{
	transcript_set ts(t, count, sid, this->single_exon_overlap);
	add(ts, mode1, mode2);
	return 0;
}

int transcript_set::add(const transcript_set &ts, int mode1, int mode2)
{
	//boost::asio::thread_pool pool(threads);
	for(auto &x : ts.mt)
	{
		map<size_t, vector<trans_item>>::iterator z = mt.find(x.first);
		if(z == mt.end())
		{
			mt.insert(x);
		}
		else
		{
			merge_sorted_trans_items(z->second, x.second, mode1, mode2, single_exon_overlap);
			/*
			vector<trans_item> &zz = z->second;
			const vector<trans_item> &xx = x.second;
			if(threads <= 0) merge_sorted_trans_items(zz, xx, mode);
			else boost::asio::post(pool, [&zz, &xx, mode] { merge_sorted_trans_items(zz, xx, mode); });
			*/
		}
	}
	//pool.join();
	return 0;
}

int transcript_set::filter(int min_count)
{
	for(auto &x: mt)
	{
		vector<trans_item> v;
		for(auto &z: x.second)
		{
			if(z.count < min_count) continue;
			v.push_back(z);
		}
		x.second = std::move(v);
	}
	return 0;
}

int transcript_set::increase_count(int count)
{
	for(auto &x : mt)
	{
		for(auto &z : x.second)
		{
			z.count += count;
		}
	}
	return 0;
}

int transcript_set::print() const
{
	printf("transcript-set: chrm = %s, mt.size() = %lu\n", chrm.c_str(), mt.size());
	return 0;
}

vector<transcript> transcript_set::get_transcripts(int min_count) const
{
	vector<transcript> v;
	for(auto &x : mt)
	{
		for(auto &z : x.second)
		{
			if(z.count < min_count) continue;
			v.push_back(z.trst);
		}
	}
	return v;
}

vector<transcript> transcript_set::get_transcripts(int min_count_single, int min_count_multi) const
{
	vector<transcript> v;
	for(auto &x : mt)
	{
		for(auto &z : x.second)
		{
			if(z.trst.exons.size() <= 1 && z.count < min_count_single) continue;
			if(z.trst.exons.size() >= 2 && z.count < min_count_multi) continue;
			v.push_back(z.trst);
		}
	}
	return v;
}

pair<bool, trans_item> transcript_set::query(const transcript &t) const
{
	pair<bool, trans_item> p;
	size_t h = t.get_intron_chain_hashing();
	auto it = mt.find(h);
	if(it == mt.end()) 
	{
		p.first = false;
		return p;
	}

	auto &v = it->second;
	for(int k = 0; k < v.size(); k++)
	{
		const transcript &x = v[k].trst;
		if(x.strand != t.strand) continue;
		bool b = x.equal1(t, single_exon_overlap);
		if(b == true) 
		{
			p.first = true;
			p.second = v[k];
			return p;
		}
	}
	p.first = true;
	return p;
}
