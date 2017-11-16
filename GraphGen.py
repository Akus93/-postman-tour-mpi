import math
from random import randint

# v1 = [1,1,1,1,2,2,2,3,4]
# v2 = [2,3,5,6,3,4,5,4,6]
# B = [0,4,8,11,14,16,18]
# v1v2 = v1 + v2
#
# v1v2.sort()
# print(v1v2)
# s = [0,2,4,6,8,9,10]
# nv1v2 = list(range(len(v1v2)))
#
#
# print(len(nv1v2))
#
# odd = 5
#
# for e in range(len(v1v2)):
#     if odd != e:
#         t = v1v2[e]
#         tm = s[t-1]
#         t2 = v1v2[e]
#         b1 = B[t2-1]
#         tm2 = ((e+1 - b1  )/2)
#         # print(tm2)
#         nv1v2[e] = math.ceil(tm + tm2)
#     else:
#         t  = v1v2[e]
#         tm = s[t -1]
#         t2 = v1v2[e]
#         b1 = B[t2-1]
#         tm2 = ((e - b1)/2)
#         # print(tm2)
#         nv1v2[e] = math.ceil( tm + tm2)
#     # print(f'dupa {e}')
# nv1v2.sort()
# print(nv1v2)
#
# nv1 = [1,1,2,2,3,4,4,6,8]
# nv2 = [3,5,9,10,6,7,9,7,10]
#
# used = 0
#
#
# #
#
# path = [0] * int(len(v1) + 1)
# odd2 = odd
# for x in range(len(v1)):
#     print(x)
#     for y in range(len(v1)):
#         if odd2 == nv1[y]:
#             path[x] = nv1[y]
#             odd2 = nv2[y]
#             nv1[y] = 0
#             nv2[y] = 0
#             break
#         elif odd2 == nv2[y]:
#             path[x] = nv2[y]
#             odd2 = nv1[y]
#             nv1[y] = 0
#             nv2[y] = 0
#             break
#
# path[len(path)-1] = odd2
# print(path)
#
import random
from datetime import datetime

def Graph_with_Euler_path(n,m):
    random.seed(datetime.now())

    nodes = []
    nodes = [n for n in range(0, n)]

    v1 = []
    v2 = []

    v1 += nodes
    v2 += nodes[1:]

    v1 += [0] * (m - n)
    v2 += [0] * (m - n + 1)
    # print(f'Wierzcholki\n {nodes}')
    v2[n - 1] = randint(0, n - 1)
    tmp = n
    tmp2 = m - n
    # print(f'V1 : \n {v1}')
    # print(f'V2 : \n {v2}')
    kv1v2 = [str(x) + ' ' + str(y) for (x, y) in zip(v1, v2)]
    v1[n - 1] = v1[0]
    v2[n - 1] = v1[0]

    tmp = m-n+1
    tmp2 = n-1
    tmp3 = [str(x) + ' ' + str(y) for (x, y) in zip(v1, v2)]

    t = v2[n-2]
    print(
        t
    )
    while (tmp):

        t2 = randint(0, n-1)
        if str(t) + ' ' + str(t2) in tmp3 or t == t2 or str(t2) + ' ' + str(t) in tmp3:
            continue
        else:
            # print(t)
            # print(t2)
            v1[tmp2] = t
            v2[tmp2] = t2
            t = t2
            tmp -= 1
            tmp2 += 1
            tmp3 = [str(x) + ' ' + str(y) for (x, y) in zip(v1, v2)]

    # print(f"mojen{n}")
    # v1[n - 1] = v1[0]
    # v2[n - 1] = v2[n - 2]
    # print(f'V1 : \n {v1}')
    # print(f'V2 : \n {v2}')
    #
    v1v2 = v1 + v2
    # print(f'V1V2:\n{v1v2}')
    v1 = [x + 1 for x in v1]
    v2 = [x + 1 for x in v2]

    kv1v2 = [str(x) + ' ' + str(y) for (x, y) in zip(v1, v2)]

    #print(kv1v2)

    file = open('graph.txt', mode='w')
    file.write("\n".join(str(item) for item in kv1v2))
    file.write("\n")


if __name__ == "__main__":
    n = 1000
    m = 20000
    if(m <= n-1 and m <=(n*(n-1)/2)):
        exit(1)
    Graph_with_Euler_path(n,m)
