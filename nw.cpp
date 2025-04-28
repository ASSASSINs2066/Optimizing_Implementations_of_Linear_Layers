#include <iostream>
#include <vector>
#include <bitset>
#include <fstream>
#include <cstdlib>
#include <pthread.h>
#include <random>
#include <sys/time.h>
#include <memory.h>

using namespace std;

#define MULTI_THREAD_FLAG 1
#define THREAD_NUM 4
#define CHOICE 33

#if CHOICE == 33
#define SIZE 32
#define FILENAME "AnS3.txt"
#else
#define SIZE 16
#endif

typedef bitset<SIZE> ROW;

typedef struct
{
    int src;
    int dst;
    bool flag;
} xpair;

typedef struct
{
    vector<xpair> seq;
    int gap;
    int start;
    int len;
} thread_data;

mt19937 rand_generator(time(NULL));

vector<ROW> get_matrix()
{
    vector<ROW> m{
        0x10821700,
        0x10100570,
        0x90080205,
        0xa0303189,
        0x00453223,
        0x20020025,
        0x831a0024,
        0x0c6090c4,
        0x0e404203,
        0x2207050c,
        0x0440c13c,
        0x00948902,
        0x8395a010,
        0x4141808a,
        0x43244820,
        0x58241085,
        0x8b044028,
        0x050604b0,
        0x028b1050,
        0x01090240,
        0xa4d0a800,
        0x58026254,
        0x04494020,
        0x8102021e,
        0xc08e01b0,
        0x00a83000,
        0x4c019909,
        0x8060ad48,
        0x23583201,
        0x16248882,
        0x000064c4,
        0x18816640

    };
    return m;
}
vector<xpair> reduce(vector<ROW> m);

bool reduce0(vector<xpair> &seq, int **table);
int reduce_step(vector<xpair> &seq);
int **get_table(vector<xpair> seq, int **table, int osize, int nsize);

int get_ones(vector<ROW> m);

vector<xpair> strgy1(vector<ROW> &m);
vector<xpair> strgy2(vector<ROW> &m);
vector<xpair> strgy3(vector<ROW> &m);
int get_ones(vector<ROW> m);
void get_trans_matrix(vector<ROW> &trans_m, vector<ROW> m);
vector<xpair> update_seq_str(vector<xpair> seq, int tab[SIZE]);

int select_oper(vector<ROW> m, vector<xpair> &max_seq, int no_reduced, int opr_type);

void build_table(vector<ROW> m, int tab[SIZE]);

vector<xpair> strgy1(vector<ROW> &m)
{
    /* compatible with non-square matrix.*/
    vector<xpair> seq;
    int row_size = m.size();
    int col_size = m[0].size();

    int *mark = new int[row_size];
    memset(mark, 0, row_size * sizeof(int));
    xpair p;
    for (unsigned col = 1; col <= col_size; ++col)
    {
        unsigned r = 0;
        while (((r < row_size) && (!m[r][col_size - col])) || (mark[r] == 1))
        {
            ++r;
        }
        if (r >= row_size)
        {
            continue;
        }
        else
        {
            mark[r] = 1;
        }

        for (unsigned i = 0; i < row_size; ++i)
        {
            if (m[i][col_size - col] && (i != r))
            {
                m[i] ^= m[r];
                p.dst = i;
                p.src = r;
                p.flag = false;
                seq.push_back(p);
            }
        }
    }
    delete[] mark;
    return seq;
}

vector<xpair> strgy2(vector<ROW> &m)
{
    vector<ROW> trans_m;
    get_trans_matrix(trans_m, m);
    vector<xpair> seq = strgy1(trans_m);
    for (int i = 0; i < seq.size(); i++)
    {
        seq[i].flag = true;
    }
    get_trans_matrix(m, trans_m);

    int tab[SIZE] = {0};
    build_table(m, tab); // get TABLE
    vector<xpair> final_seq = update_seq_str(seq, tab);
    return final_seq;
}
vector<xpair> strgy3(vector<ROW> &m)
{
    vector<xpair> tmp_seq;
    while (get_ones(m) != m.size())
    {
        vector<xpair> base_oper;
        base_oper.clear();
        int value = 0;

        bool changed = false; // 标志变量，检测矩阵是否发生变化

        int Reduced_Row = select_oper(m, base_oper, value, 0);

        vector<ROW> trans_m;
        get_trans_matrix(trans_m, m);

        int Reduced_Col = select_oper(trans_m, base_oper, Reduced_Row, 1);
        if (base_oper.size() >= 1)
        {
            int rand_num = rand_generator() % base_oper.size();
            // cout<<rand_num<<endl;
            if (!base_oper[rand_num].flag)
            {
                ROW old_row = m[base_oper[rand_num].dst]; // 保存原始行
                m[base_oper[rand_num].dst] ^= m[base_oper[rand_num].src];
                if (m[base_oper[rand_num].dst] != old_row) // 检查是否发生变化
                {
                    changed = true;
                }
                tmp_seq.push_back(base_oper[rand_num]);
            }
            else if (base_oper[rand_num].flag)
            {
                vector<ROW> trans_tmp_m;
                get_trans_matrix(trans_tmp_m, m);
                ROW old_row = trans_tmp_m[base_oper[rand_num].dst]; // 保存原始行
                trans_tmp_m[base_oper[rand_num].dst] ^= trans_tmp_m[base_oper[rand_num].src];
                if (trans_tmp_m[base_oper[rand_num].dst] != old_row) // 检查是否发生变化
                {
                    changed = true;
                }
                get_trans_matrix(m, trans_tmp_m);
                tmp_seq.push_back(base_oper[rand_num]);
            }
        }
        else
        {
            break;
        }

        if (!changed) // 如果矩阵没有发生变化，退出循环
        {
            break;
        }
    }

    if (get_ones(m) != m.size())
    {
        int rnd = rand_generator() % 2;
        if (rnd == 0)
        {
            vector<xpair> seq_2(strgy1(m));
            tmp_seq.insert(tmp_seq.end(), seq_2.begin(), seq_2.end());
        }
        else if (rnd == 1)
        {
            vector<xpair> seq_2(strgy2(m));
            tmp_seq.insert(tmp_seq.end(), seq_2.begin(), seq_2.end());
        }
    }

    int tab[SIZE] = {0};
    build_table(m, tab);
    vector<xpair> final_seq = update_seq_str(tmp_seq, tab);
    return final_seq;
}
// vector<xpair> strgy3(vector<ROW> &m)
// {
//     vector<xpair> tmp_seq;
//     while (get_ones(m) != m.size())
//     {
//         vector<xpair> base_oper;
//         base_oper.clear();
//         int value = 0;

//         int Reduced_Row = select_oper(m, base_oper, value, 0);

//         vector<ROW> trans_m;
//         get_trans_matrix(trans_m, m);

//         int Reduced_Col = select_oper(trans_m, base_oper, Reduced_Row, 1);
//         if (base_oper.size() >= 1)
//         {
//             int rand_num = rand_generator() % base_oper.size();

//             if (!base_oper[rand_num].flag)
//             {
//                 m[base_oper[rand_num].dst] ^= m[base_oper[rand_num].src];
//                 tmp_seq.push_back(base_oper[rand_num]);
//             }

//             else if (base_oper[rand_num].flag)
//             {
//                 vector<ROW> trans_tmp_m;
//                 get_trans_matrix(trans_tmp_m, m);
//                 trans_tmp_m[base_oper[rand_num].dst] ^= trans_tmp_m[base_oper[rand_num].src];
//                 get_trans_matrix(m, trans_tmp_m);
//                 tmp_seq.push_back(base_oper[rand_num]);
//             }
//         }
//         else
//         {
//             break;
//         }
//     }

//     if (get_ones(m) != m.size())
//     {
//         // int rnd = rand()%2;
//         int rnd = rand_generator() % 2;
//         if (rnd == 0)
//         {
//             vector<xpair> seq_2(strgy1(m));
//             tmp_seq.insert(tmp_seq.end(), seq_2.begin(), seq_2.end());
//         }
//         else if (rnd == 1)
//         {
//             vector<xpair> seq_2(strgy2(m));
//             tmp_seq.insert(tmp_seq.end(), seq_2.begin(), seq_2.end());
//         }
//     }

//     int tab[SIZE] = {0};
//     build_table(m, tab);
//     vector<xpair> final_seq = update_seq_str(tmp_seq, tab);
//     return final_seq;
// }

int select_oper(vector<ROW> m, vector<xpair> &max_seq, int no_reduced, int opr_type)
{
    int no_before = 0;
    int no_after = 0;
    xpair new_ele;

    for (int i = 0; i < m.size(); i++)
    {
        vector<ROW> tmp_m(m);
        for (int j = 0; j < m.size(); j++)
        {
            if (i != j)
            {
                no_before = tmp_m[j].count();
                tmp_m[j] ^= tmp_m[i];
                no_after = tmp_m[j].count();
                if ((no_before - no_after) > 0)
                {
                    if (no_reduced < (no_before - no_after))
                    {
                        no_reduced = no_before - no_after;
                        max_seq.clear();
                        new_ele.src = i;
                        new_ele.dst = j;
                        if (opr_type == 0)
                        {
                            new_ele.flag = false;
                        }
                        else if (opr_type == 1)
                        {
                            new_ele.flag = true;
                        }
                        max_seq.push_back(new_ele);
                    }

                    else if (no_reduced == (no_before - no_after))
                    {
                        new_ele.src = i;
                        new_ele.dst = j;
                        if (opr_type == 0)
                        {
                            new_ele.flag = false;
                        }
                        else if (opr_type == 1)
                        {
                            new_ele.flag = true;
                        }
                        max_seq.push_back(new_ele);
                    }
                }
            }
        }
    }
    return no_reduced;
}

int get_ones(vector<ROW> m)
{
    int s = 0;
    for (int i = 0; i < m.size(); i++)
        s += m[i].count();
    return s;
}

void get_trans_matrix(vector<ROW> &trans_m, vector<ROW> m)
{
    trans_m.clear();
    for (int i = 0; i < m.size(); i++)
    {
        bitset<SIZE> tmp(0);
        for (int j = 0; j < m.size(); j++)
        {
            tmp[m[0].size() - 1 - j] = m[j][m[0].size() - i - 1];
        }
        trans_m.push_back(tmp);
    }
}

void build_table(vector<ROW> m, int tab[SIZE])
{
    int ind = 0;
    for (int i = 0; i < m.size(); i++)
    {
        for (int j = 0; j < m.size(); j++)
        {
            if (m[i].test(j))
            {
                ind = m.size() - 1 - j;
                break;
            }
        }
        tab[ind] = i;
    }
}

vector<xpair> update_seq_str(vector<xpair> seq, int tab[SIZE])
{
    vector<xpair> tmp_seq;
    for (int i = 0; i < seq.size(); i++)
    {
        if (!seq[i].flag)
        {
            tmp_seq.push_back(seq[i]);
        }
    }

    for (int i = seq.size() - 1; i >= 0; i--)
    {

        if (seq[i].flag)
        {
            xpair new_ele;
            new_ele.src = tab[seq[i].dst];
            new_ele.dst = tab[seq[i].src];
            new_ele.flag = false;
            tmp_seq.push_back(new_ele);
        }
    }
    return tmp_seq;
}

vector<ROW> get_reduced_matrix(vector<xpair> seq, vector<ROW> m)
{
    vector<ROW> tmp_m(m);
    for (int i = 0; i < seq.size(); i++)
        tmp_m[seq[i].dst] ^= tmp_m[seq[i].src];
    return tmp_m;
}

vector<ROW> get_identity()
{
    vector<ROW> m;
    ROW tmp(0);
    for (int i = 0; i < SIZE; i++)
    {
        m.push_back(tmp);
        m[i][SIZE - 1 - i] = 1;
    }
    return m;
}

bool exchange(xpair a, xpair b)
{
    if (a.dst == b.dst)
        return true;
    if (a.src == b.src)
        return true;
    if ((a.src != b.dst) && (a.dst != b.src))
        return true;
    return false;
}

int **get_table(vector<xpair> seq, int **table, int osize, int nsize)
{
    if (table != NULL)
    {
        for (int i = 0; i < osize; i++)
            delete[] table[i];
        delete[] table;
    }
    table = new int *[nsize];
    for (int i = 0; i < nsize; i++)
        table[i] = new int[nsize];

    int s = seq.size();
    for (int i = 0; i < s; i++)
    {
        table[i][i] = 1;
        for (int j = i + 1; j < s; j++)
        {
            if (exchange(seq[i], seq[j]))
                table[i][j] = 1;
            else
            {
                for (int k = j; k < s; k++)
                    table[i][k] = 0;
                break;
            }
        }
        for (int j = i - 1; j >= 0; j--)
        {
            if (exchange(seq[i], seq[j]))
                table[i][j] = 1;
            else
            {
                for (int k = j; k >= 0; k--)
                    table[i][k] = 0;
                break;
            }
        }
    }
    return table;
}

void delete_table(int **table, int size)
{
    if (table != NULL)
    {
        for (int i = 0; i < size; i++)
            delete[] table[i];
        delete[] table;
    }
}

bool exchange_set(vector<xpair> seq, int start, int end, xpair p)
{
    if (start <= end)
    {
        for (int i = start; i <= end; i++)
        {
            if (!exchange(p, seq[i]))
                return false;
        }
    }
    else
    {
        for (int i = start; i >= end; i--)
        {
            if (!exchange(p, seq[i]))
                return false;
        }
    }
    return true;
}

bool reduce0(vector<xpair> &seq, int **table)
{
    int s = seq.size();
    int index;

    for (int i = 0; i < s; i++)
    {
        for (int j = i + 1; j < s; j++)
        {
            if ((table[i][j - 1] == 0) && (table[j][i + 1] == 0))
                break;
            if (seq[i].src != seq[j].src)
                continue;
            for (int k = j + 1; k < s; k++)
            {
                if ((table[k][j + 1] == 0) && (table[j][k - 1] == 0))
                    break;
                if (seq[k].dst != seq[i].dst)
                {
                    if (seq[k].dst != seq[j].dst)
                        continue;
                    else
                    {
                        if (seq[k].src != seq[i].dst)
                            continue;
                    }
                }
                else
                {
                    if (seq[k].src != seq[j].dst)
                        continue;
                }
                if ((table[i][j - 1]) && table[k][j + 1])
                    index = j - 1; // fixpoint = j;
                else if (table[j][i + 1] && table[k][j + 1] && exchange_set(seq, j - 1, i + 1, seq[k]))
                    index = i; // fixpoint = i;
                else if (table[i][j - 1] && table[j][k - 1] && exchange_set(seq, j + 1, k - 1, seq[i]))
                    index = k - 2; // fixpoint = k;
                else
                    continue;
                xpair tmp1, tmp2;
                tmp1 = seq[k];
                tmp2.src = seq[i].src;
                tmp2.dst = seq[k].src;
                seq.erase(seq.begin() + k);
                seq.erase(seq.begin() + j);
                seq.erase(seq.begin() + i);
                seq.insert(seq.begin() + index, tmp1);
                seq.insert(seq.begin() + index + 1, tmp2);
                return true;
            }
        }
    }
    return false;
}

bool reduce1(vector<xpair> &seq, int **table)
{
    int s = seq.size();
    int index;

    for (int i = 0; i < s; i++)
    {
        for (int j = i + 1; j < s; j++)
        {
            if ((table[i][j - 1] == 0) && (table[j][i + 1] == 0))
                break;
            if ((seq[i].src != seq[j].dst) && (seq[i].dst != seq[j].dst))
                continue;
            for (int k = j + 1; k < s; k++)
            {
                if ((table[k][j + 1] == 0) && (table[j][k - 1] == 0))
                    break;
                if (seq[k].src != seq[j].src)
                    continue;
                else
                {
                    if (seq[j].dst == seq[i].dst)
                    {
                        if (seq[k].dst != seq[i].src)
                            continue;
                    }
                    if (seq[j].dst == seq[i].src)
                    {
                        if (seq[k].dst != seq[i].dst)
                            continue;
                    }
                }
                if ((table[i][j - 1]) && table[k][j + 1])
                    index = j - 1; // fixpoint = j;
                else if (table[j][i + 1] && table[k][j + 1] && exchange_set(seq, j - 1, i + 1, seq[k]))
                    index = i; // fixpoint = i;
                else if (table[i][j - 1] && table[j][k - 1] && exchange_set(seq, j + 1, k - 1, seq[i]))
                    index = k - 2; // fixpoint = k;
                else
                    continue;
                xpair tmp1, tmp2;
                tmp1.src = seq[j].src;
                tmp1.dst = seq[i].src;
                tmp2 = seq[i];
                seq.erase(seq.begin() + k);
                seq.erase(seq.begin() + j);
                seq.erase(seq.begin() + i);
                seq.insert(seq.begin() + index, tmp1);
                seq.insert(seq.begin() + index + 1, tmp2);
                return true;
            }
        }
    }
    return false;
}

bool reduce2(vector<xpair> &seq, int **table)
{
    int s = seq.size();
    int index;

    for (int i = 0; i < s; i++)
    {
        for (int j = i + 1; j < s; j++)
        {
            if ((table[i][j - 1] == 0) && (table[j][i + 1] == 0))
                break;
            if ((seq[j].src != seq[i].dst) && (seq[j].dst != seq[i].dst))
                continue;
            for (int k = j + 1; k < s; k++)
            {
                if ((table[k][j + 1] == 0) && (table[j][k - 1] == 0))
                    break;
                if (seq[k].src != seq[i].src)
                    continue;
                else
                {
                    if (seq[j].src == seq[i].dst)
                    {
                        if (seq[k].dst != seq[j].dst)
                            continue;
                    }
                    if (seq[j].dst == seq[i].dst)
                    {
                        if (seq[k].dst != seq[j].src)
                            continue;
                    }
                }
                if ((table[i][j - 1]) && table[k][j + 1])
                    index = j - 1; // fixpoint = j;
                else if (table[j][i + 1] && table[k][j + 1] && exchange_set(seq, j - 1, i + 1, seq[k]))
                    index = i; // fixpoint = i;
                else if (table[i][j - 1] && table[j][k - 1] && exchange_set(seq, j + 1, k - 1, seq[i]))
                    index = k - 2; // fixpoint = k;
                else
                    continue;
                xpair tmp1, tmp2;
                tmp1.dst = seq[k].dst;
                tmp2.dst = seq[i].dst;
                if (seq[j].src == seq[i].dst)
                {
                    tmp1.src = seq[j].src;
                    tmp2.src = seq[i].src;
                }
                else
                {
                    tmp1.src = seq[i].src;
                    tmp2.src = seq[j].src;
                }
                seq.erase(seq.begin() + k);
                seq.erase(seq.begin() + j);
                seq.erase(seq.begin() + i);
                seq.insert(seq.begin() + index, tmp1);
                seq.insert(seq.begin() + index + 1, tmp2);
                return true;
            }
        }
    }
    return false;
}

bool reduce3(vector<xpair> &seq, int **table)
{
    int s = seq.size();
    int index;

    for (int i = 0; i < s; i++)
    {
        for (int j = i + 1; j < s; j++)
        {
            if ((table[i][j - 1] == 0) && (table[j][i + 1] == 0))
                break;
            if (seq[j].dst != seq[i].src)
                continue;
            if (seq[j].src != seq[i].dst)
                continue;
            if (table[i][j - 1])
                index = j - 1;
            else
                index = i;
            xpair tmp = seq[j];
            seq.erase(seq.begin() + j);
            seq.erase(seq.begin() + i);
            seq.insert(seq.begin() + index, tmp);
            int mask = seq[index].dst ^ seq[index].src;
            for (int k = index + 1; k < s; k++)
            {
                if ((seq[k].dst == seq[index].dst) || (seq[k].dst == seq[index].src))
                    seq[k].dst ^= mask;
                if ((seq[k].src == seq[index].dst) || (seq[k].src == seq[index].src))
                    seq[k].src ^= mask;
            }
            return true;
        }
    }
    return false;
}

int reduce_step(vector<xpair> &seq)
{
    int **tab = NULL;
    tab = get_table(seq, tab, 0, seq.size());
    int i = 0;
    int NUM = 4;
    int counter = 0;
    while (counter != NUM)
    {
        switch (i)
        {
        case 0:
        {
            if (reduce0(seq, tab))
            {
                // cout << "activate case 0"<<endl;
                tab = get_table(seq, tab, seq.size() + 1, seq.size());
                // i = 0;
                counter = 0;
            }
            else
            {
                i = (i + 1) % NUM;
                counter++;
            }
            break;
        }
        case 1:
        {
            if (reduce1(seq, tab))
            {
                // cout << "activate case 1" << endl;
                tab = get_table(seq, tab, seq.size() + 1, seq.size());
                // i = 1;
                counter = 0;
            }
            else
            {
                i = (i + 1) % NUM;
                counter++;
            }
            break;
        }
        case 2:
        {
            if (reduce2(seq, tab))
            {
                // cout << "activate case 2" <<endl;
                tab = get_table(seq, tab, seq.size() + 1, seq.size());
                counter = 0;
            }
            else
            {
                i = (i + 1) % NUM;
                counter++;
            }
            break;
        }
        case 3:
        {

            if (reduce3(seq, tab))
            {
                // cout << "activate case 3" << endl;
                tab = get_table(seq, tab, seq.size() + 1, seq.size());
                counter = 0;
            }
            else
            {
                // cout << "case 3" << endl;
                i = (i + 1) % NUM;
                counter++;
            }
        }
        }
    }
    delete_table(tab, seq.size());
    return seq.size();
}

void update_seq(vector<xpair> &seq, int tab[SIZE], int start)
{

    for (int i = start; i < seq.size(); i++)
    {
        seq[i].src = tab[seq[i].src];
        seq[i].dst = tab[seq[i].dst];
    }
}

void get_equivalent_seq(vector<xpair> &seq, int gap, int start)
{

    vector<ROW> m;
    bitset<SIZE> tmp;
    for (int i = 0; i < SIZE; i++)
    {
        tmp = 0;
        tmp[SIZE - 1 - i] = 1;
        m.push_back(tmp);
    }
    for (int i = start + gap - 1; i >= start; i--)
    {
        m[seq[i].dst] ^= m[seq[i].src];
    }

    vector<xpair> seq_here(strgy3(m));
    seq.erase(seq.begin() + start, seq.begin() + start + gap);
    seq.insert(seq.begin() + start, seq_here.begin(), seq_here.end());
    int tab[SIZE] = {0};
    build_table(m, tab);
    update_seq(seq, tab, start + seq_here.size());
}

void *reduce_thread(void *d)
{
    thread_data *data = (thread_data *)d;
    get_equivalent_seq(data->seq, data->gap, data->start);
    data->len = reduce_step(data->seq);
    pthread_exit(NULL);
}

bool check_consistence(vector<xpair> seq)
{
    vector<ROW> m = get_matrix();
    for (int i = 0; i < seq.size(); i++)
        m[seq[i].dst] ^= m[seq[i].src];
    if (get_ones(m) == SIZE)
        return true;
    else
        return false;
}

vector<xpair> reduce(vector<ROW> m)
{

    vector<ROW> tmp_m(m);

    vector<xpair> seq = strgy3(m);
    int l = reduce_step(seq);

    int gap = l;
    int start = 0;

    thread_data data[THREAD_NUM];
    bool flag = true;
    while (flag)
    {
        for (int j = 0; j < THREAD_NUM; j++)
        {
            data[j].seq = seq;
            data[j].gap = gap;
            data[j].start = start;
            if ((gap != 3) || (start != l - gap))
            {
                if (start == l - gap)
                {
                    gap--;
                    start = 0;
                }
                else
                {
                    start++;
                }
            }
            else
            {
                flag = false;
            }
        }
        int rc;
        pthread_t threads[THREAD_NUM];
        pthread_attr_t attr;
        void *status;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        for (int j = 0; j < THREAD_NUM; j++)
        {

            rc = pthread_create(&threads[j], NULL, reduce_thread, (void *)(data + j));
            if (rc)
            {
                cout << "Error: unable to create thread!" << endl;
                exit(-1);
            }
        }
        // pthread_attr_destroy(&attr);
        for (int j = 0; j < THREAD_NUM; j++)
        {
            rc = pthread_join(threads[j], &status);
            if (rc)
            {
                cout << "Error: unable to join, " << rc << endl;
                exit(-1);
            }
        }
        for (int j = 0; j < THREAD_NUM; j++)
        {
            if (data[j].len < l)
            {
                seq = data[j].seq;
                l = data[j].len;
                gap = l;
                start = 0;
            }
        }
        // cout << "gap = " << gap << endl;
    }
    return seq;
}

void print_seq(vector<ROW> m, vector<ROW> tmp_m, vector<xpair> seq);
int main()
{
    int counter = 100000;
    struct timeval start;
    struct timeval end;
    gettimeofday(&start, NULL);

    int num = 0;
    cout << FILENAME << endl;
    // while((num++) < 1)
    // srand(time(0));
    while (true)
    {
        vector<ROW> m = get_matrix();

        vector<xpair> seq = reduce(m);

        vector<ROW> tmp_m = get_reduced_matrix(seq, m);

        cout << "minimal = " << counter << "   current = " << seq.size() << endl
             << endl;

        if (seq.size() < counter)
        {
            counter = seq.size();
            print_seq(m, tmp_m, seq);
        }
    }

    gettimeofday(&end, NULL);
    float time = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    time /= 1000000;
    cout << "Time used = " << time << endl;

    return 0;
}

void print_seq(vector<ROW> m, vector<ROW> tmp_m, vector<xpair> seq)
{
    ofstream f;
    f.open(FILENAME);
    f << "Original Matrix:" << endl;
    for (int i = 0; i < m.size(); i++)
        f << m[i] << endl;
    f << endl
      << endl;
    f << "Reduced Matrix:" << endl;
    for (int i = 0; i < tmp_m.size(); i++)
        f << tmp_m[i] << endl;
    f << endl
      << endl;

    f << "Xor Count = " << seq.size() << endl;

    int tab[SIZE] = {0};
    for (int i = 0; i < tmp_m.size(); i++)
    {
        for (int j = 0; j < tmp_m.size(); j++)
        {
            if (tmp_m[i].test(tmp_m.size() - 1 - j))
            {
                tab[i] = j;
                break;
            }
        }
    }

    for (int i = seq.size() - 1; i >= 0; i--)
    { // f << seq[i].dst << " <- " << seq[i].src << endl;
        f << "x[" << tab[seq[i].dst] << "] = x[" << tab[seq[i].dst] << "] ^ x[" << tab[seq[i].src] << "]";
        tmp_m[seq[i].dst] ^= tmp_m[seq[i].src];
        bool flag = false;
        for (int j = 0; j < m.size(); j++)
        {
            if (tmp_m[seq[i].dst] == m[j])
            {
                flag = true;
                f << "    y[" << j << "]" << endl;
                break;
            }
        }

        if (!flag)
        {
            f << endl;
        }
    }
    f.close();
}