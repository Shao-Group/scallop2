/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop2
(c) 2021 by  Qimin Zhang, Mingfu Shao, and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "fragment.h"
#include "util.h"
#include "config.h"
#include <cstdio>

fragment::fragment(hit *x1, hit *x2)
	: h1(x1), h2(x2)
{
	frag_type = 1;
	is_compatible = 0;
	pi = -1;
	fidx = -1;
	fake_hit_index = -1;
	HS_frag = false;
	b1 = false;
	b2 = false;
	k1l = 0;
	k1r = 0;
	k2l = 0;
	k2r = 0;
	lpos = -1;
	rpos = -1;
	cnt = 1;
	paths.clear();
}

int fragment::print(int index)
{
	printf("fragment %d: name = %s, cnt = %d, lpos = %d, rpos = %d, len = %d, k1l = %d, k1r = %d, k2l = %d, k2r = %d, b1 = %c, b2 = %c, v1 = %lu, v2 = %lu, #paths = %lu, frag_type = %d, pi = %d, fidx = %d,type = %d, is_compatible = %d, fake_hit_index = %d\n",
			index, h1->qname.c_str(), cnt, lpos, rpos, rpos - lpos, k1l, k1r, k2l, k2r, b1 ? 'T' : 'F', b2 ? 'T' : 'F', h1->vlist.size(), h2->vlist.size(), paths.size(), frag_type, pi, fidx, type, is_compatible, fake_hit_index);
	h1->print();
	h2->print();

	printf(" v1 = ( ");
	printv(decode_vlist(h1->vlist));
	printf(") v2 = ( ");
	printv(decode_vlist(h2->vlist));
	printf(")\n");

	for(int k = 0; k < paths.size(); k++)
	{
		printf("\t\t");
		paths[k].print_bridge(k);
	}
	return 0;
}

int fragment::clear()
{
	frag_type = 1;
	is_compatible = 0;
	b1 = false;
	b2 = false;
	k1l = 0;
	k1r = 0;
	k2l = 0;
	k2r = 0;
	lpos = -1;
	rpos = -1;
	cnt = 1;
	return 0;
}

int fragment::append(const fragment &f) 
{
	assert(f.cnt == 1);
	hit *x1 = f.h1;
	hit *x2 = f.h2;
	assert(x1->next == NULL);
	assert(x2->next == NULL);
	x1->next = h1;
	x2->next = h2;
	h1 = x1;
	h2 = x2;
	cnt++;
	return 0;
}

bool fragment::equal(const fragment &f) const
{
	if(h1->vlist != f.h1->vlist) return false;
	if(h2->vlist != f.h2->vlist) return false;
	if(lpos != f.lpos) return false;
	if(rpos != f.rpos) return false;
	if(k1l != f.k1l) return false;
	if(k1r != f.k1r) return false;
	if(k2l != f.k2l) return false;
	if(k2r != f.k2r) return false;
	if(b1 != f.b1) return false;
	if(b2 != f.b2) return false;
	return true;
}

int fragment::set_paired(bool b)
{
	hit *x = h1;
	hit *y = h2;
	while(x != NULL && y != NULL)
	{
		x->paired = b;
		y->paired = b;
		x = x->next;
		y = y->next;
	}
	return 0;
}

int fragment::set_bridged(bool b)
{
	hit *x = h1;
	hit *y = h2;
	while(x != NULL && y != NULL)
	{
		//if(x->bridged == false) x->bridged = b;
		//if(y->bridged == false) y->bridged = b;
		x->bridged = b;
		y->bridged = b;
		x = x->next;
		y = y->next;
	}
	return 0;
}

bool compare_fragment(const fragment &f1, const fragment &f2)
{
	if(f1.h1->vlist.size() < f2.h1->vlist.size()) return true;
	if(f1.h1->vlist.size() > f2.h1->vlist.size()) return false;
	if(f1.h2->vlist.size() < f2.h2->vlist.size()) return true;
	if(f1.h2->vlist.size() > f2.h2->vlist.size()) return false;

	for(int k = 0; k < f1.h1->vlist.size(); k++)
	{
		if(f1.h1->vlist[k] < f2.h1->vlist[k]) return true;
		if(f1.h1->vlist[k] > f2.h1->vlist[k]) return false;
	}

	for(int k = 0; k < f1.h2->vlist.size(); k++)
	{
		if(f1.h2->vlist[k] < f2.h2->vlist[k]) return true;
		if(f1.h2->vlist[k] > f2.h2->vlist[k]) return false;
	}

	if(f1.b1 == true && f2.b1 == false) return true;
	if(f1.b1 == false && f2.b1 == true) return false;
	if(f1.b2 == true && f2.b2 == false) return true;
	if(f1.b2 == false && f2.b2 == true) return false;

	if(f1.k1l + f1.k1r + f1.k2l + f1.k2r < f2.k1l + f2.k1r + f2.k2l + f2.k2r) return true;
	if(f1.k1l + f1.k1r + f1.k2l + f1.k2r > f2.k1l + f2.k1r + f2.k2l + f2.k2r) return false;
	if(f1.k2l + f1.k1r < f2.k2l + f2.k1r) return true;
	if(f1.k2l + f1.k1r > f2.k2l + f2.k1r) return false;
	if(f1.k1l + f1.k2r < f2.k1l + f2.k2r) return true;
	if(f1.k1l + f1.k2r > f2.k1l + f2.k2r) return false;

	return (f1.lpos < f2.lpos);
}
