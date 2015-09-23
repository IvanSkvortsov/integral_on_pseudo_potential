#ifndef __ECP_PRIMITIVE_H__
#define __ECP_PRIMITIVE_H__

template<class T>
struct ecp_primitive
{
	T alp, r[3];
	int n, l;
	ecp_primitive():alp(0), n(0), l(0)
	{
		for(int i = 0; i < 3; ++i) r[i] = T(0);
	}
};

#endif//__ECP_PRIMITIVE_H__
