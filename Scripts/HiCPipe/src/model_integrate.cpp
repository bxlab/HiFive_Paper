#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <stdarg.h>
#include <stdlib.h>
#include <algorithm>

// Access matrix in vector form
#define MATRIX_INDEX(x, y, dim) ((x-1) + (y-1)*dim)

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/*

 Input
 -----
 1. A table of fragment ends (fends). Each fend has a list of features
    (such as frag_len_bin, gc_bin, etc.)

 2. A list of 0 or more model features. Each model feature is defined by:
    - field name: same name as the one in the fend table
    - file: contains correction matrix for feature
    - size: matrix dimension

 3. Two list (for dimension x and y) of 1 or more output features. Each feature is defined by:
    - field name: same name as the one in the fend table
    - from: feature start index
    - to: feature to index

 4. A prior on the interaction probablity

 Output
 ------
 The program outputs a single file that contains expected counts (according to model)
 for the requested output bins.

 Examples
 --------
 The following counts the size of the bins:
 %> expected_count fends.table 1 0 1 coord_bin 1 1000 1 coord_bin 1 1000

 The following checks model performance (fragment length only) in gc plane:
 %> expected_count fends.table 1 1 frag_len_bin frag_len.f 50 1 gc_bin 1 50 1 gc_bin 1 50

 Remarks
 -------
 All input and output indices are 1-based. The model matrix is represented as a vector,
 so it is zero based.
 
*/
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures
//////////////////////////////////////////////////////////////////////////////////////////////////

// The model has a set of features. Each feature is factor which is a function of
// the two interacting fends 
struct ModelFeatures
{
    int size;
    vector<string> fields;
    vector<string> fns;
    vector<int> sizes;
};

// We integrate over two sets of bins, one for each dimension (x,y)
struct OutputFeatures
{
    int size;
    vector<string> fields;
    vector<int> from;
    vector<int> to;
};

// every fend is defined by a set of levels, one for each model feature
struct Fend
{
    int index;
    int chrom_code;
    int coord;
    vector<int> levels;
};

// every bin holds a collection of all fends that have the same levels for all model features
struct Bin
{
    vector<Fend> fends;
};

// filter on contact type
enum Filter { f_both, f_trans, f_close_cis, f_far_cis };

//////////////////////////////////////////////////////////////////////////////////////////////////
// Low level utility functions
//////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
string join_strings(T it, T end, string delim, string suffix)
{
    string result = *it + suffix;
    
    while (++it != end)
        result += delim + (*it) + suffix;
    return result;
}

template <typename T>
string NumberToString ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}

void massert(bool cond, char *fmt, ...)
{
    if (cond)
        return;
    
    fprintf(stderr, "Error: ");
    
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr, fmt, argp);
    va_end(argp);
    
    fprintf(stderr, "\n");
    exit(-1);
}

// use a global structure for chromosome codes
static map<string, int> chrom_codes;
static map<int, string> rev_chrom_codes;
static int next_code = 1;
int chr2int(string chr)
{
    if (chrom_codes.find(chr) == chrom_codes.end())
    {
        chrom_codes[chr] = next_code;
        rev_chrom_codes[next_code] = chr;
        next_code++;
    }
    return chrom_codes[chr];
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// I/O utility functions
//////////////////////////////////////////////////////////////////////////////////////////////////

void split_line(istream &in, vector<string> &fields, char delim)
{
	fields.resize(0);
	string field;
	while(in) {
		char c = in.get();
		if(c == '\r') {
			continue;
		}
		if(c == '\n') {
			fields.push_back(field);
			break;
		}
		if(c == delim) {
			fields.push_back(field);
			field.resize(0);
		} else {
			field.push_back(c);
		}
	}
}

string levels_to_string(const vector<int>& levels)
{
    assert(levels.size() > 0);
    string result = NumberToString(levels[0]);
    for (unsigned int i = 1; i<levels.size(); i++)
        result += string("\t") + NumberToString(levels[i]);
    return result;
}

void parse_levels(vector<int>& levels, const vector<int>& fields_ind, const vector<string>& fields)
{
    levels.resize(fields_ind.size());
    for (unsigned int i = 0; i<fields_ind.size(); i++)
        levels[i] = atoi(fields[fields_ind[i]].c_str());
}

void add_fend_to_obins(const Fend& fend,
                       map< string, Bin > &obins,
                       const vector<int>& obin_fields_ind,
                       const vector<string>& fields)
{
    vector<int> obin_levels;
    parse_levels(obin_levels, obin_fields_ind, fields);
    vector<Fend>& fends = obins[levels_to_string(obin_levels)].fends;
    fends.push_back(fend);
}

void init_field_indices(vector<int>& field_ind, const vector<string> &fields, const vector<string>& titles)
{
    field_ind.resize(fields.size(), -1);
    for (unsigned int i = 0; i<titles.size(); i++)
        for (unsigned int j = 0; j<fields.size(); j++)
            if (titles[i] == fields[j])
                field_ind[j] = i;
    for (unsigned int i = 0; i<field_ind.size(); i++)
        massert(field_ind[i] != -1, "field %s not found", fields[i].c_str());
}

int get_field_index(string field, const vector<string>& titles)
{
    int result = -1;
    for (unsigned int i = 0; i<titles.size(); i++)
        if (titles[i] == field)
            result = i;
    assert(result != -1);
    return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Read fend file
//////////////////////////////////////////////////////////////////////////////////////////////////

void read_fends_table(const string& fn, const ModelFeatures& m_features,
                      const OutputFeatures& o_features_x, const OutputFeatures& o_features_y,
                      const int from_fend, const int to_fend, 
                      map<string, Bin>& obins_x, map<string, Bin>& obins_y)
{
    ifstream in(fn.c_str());
    massert(in.is_open(), "could not open file %s", fn.c_str());
	vector<string> fields;
    char delim = '\t';

    vector<int> m_features_ind(m_features.size, -1);
    vector<int> o_features_x_ind(o_features_x.size, -1);
    vector<int> o_features_y_ind(o_features_y.size, -1);

    // parse header
    split_line(in, fields, delim);
    init_field_indices(m_features_ind, m_features.fields, fields);
    init_field_indices(o_features_x_ind, o_features_x.fields, fields);
    init_field_indices(o_features_y_ind, o_features_y.fields, fields);

    int fend_ind = get_field_index("fend", fields);
    int chr_ind = get_field_index("chr", fields);
    int coord_ind = get_field_index("coord", fields);
    
	while(in) {
		split_line(in, fields, delim);
		if(fields.size() == 0)
			return;
        
        Fend fend;
        fend.index = atoi(fields[fend_ind].c_str());
        fend.chrom_code = chr2int(fields[chr_ind]);
        fend.coord = atoi(fields[coord_ind].c_str());
        
        // classify fend according to its model levels
        parse_levels(fend.levels, m_features_ind, fields);

        // add fend to output bins according to the output levels
        add_fend_to_obins(fend, obins_x, o_features_x_ind, fields);

        if ( (from_fend == 0 && to_fend == 0) || ( fend.index >= from_fend && fend.index <= to_fend) )
            add_fend_to_obins(fend, obins_y, o_features_y_ind, fields);
	}
    
    in.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Read feature file
//////////////////////////////////////////////////////////////////////////////////////////////////

void read_feature_table(string fn, vector<double>& result, int mat_dim)
{
    ifstream in(fn.c_str());
    massert(in.is_open(), "could not open file %s", fn.c_str());
	vector<string> fields;
    char delim = '\t';
    const unsigned int width = 3;
    bool first = true;
	while(in) {
		split_line(in, fields, delim);
		if(fields.size() == 0)
			return;
        if (first)
        {
            first = false;
            continue;
        }
        
		assert(fields.size() == width);
        
        int i = atoi(fields[0].c_str());
        int j = atoi(fields[1].c_str());
        double prob = atof(fields[2].c_str());
        
        int v1 = MATRIX_INDEX(i, j, mat_dim);
        int v2 = MATRIX_INDEX(j, i, mat_dim);

        result[v1] = prob;
        result[v2] = prob;
	}
    
    in.close();
}

void read_feature_tables(const vector<string>& feature_fns, const vector<int>& sizes,
                      vector < vector<double> >& features)
{
    for (unsigned int i=0; i < feature_fns.size(); i++)
    {
        vector<double>& feature = features[i];
        feature.resize(sizes[i]*sizes[i]);
        read_feature_table(feature_fns[i], feature, sizes[i]);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Debug
//////////////////////////////////////////////////////////////////////////////////////////////////

void dump_bins(map<string, Bin>& bins, string title)
{    
    cerr << title << " bins" << endl;
    map<string, Bin>::iterator it = bins.begin();
    while(it != bins.end())
    {
        cerr << "Bin: " << (*it).first << endl;

        vector<Fend>& fends = (*it).second.fends;
        vector<Fend>::iterator jt = fends.begin();
        while(jt != fends.end())
        {
            Fend& fend = *jt;
            cerr << " " << fend.index << " " << fend.chrom_code << " " << rev_chrom_codes[fend.chrom_code] << ": ";
            for (unsigned int i = 0; i < fend.levels.size(); i++)
                cerr << fend.levels[i] << " ";
            cerr << endl;
            jt++;
        }
        it++;
        cerr << "=============================================================" << endl;
    }
}

void dump_feature(const vector<double>& f, int mat_dim)
{
    cerr << "==============================================" << endl;
    cerr << "i\tj\tprob" << endl;
    for (unsigned int k = 0; k < f.size(); k++)
    {
        
        int i = k / mat_dim + 1;
        int j = k - (i-1) * mat_dim + 1;
        
        double prob = f[k];
        cerr << i << "\t" << j << "\t" << prob << endl;
    }
}

void dump_features(const vector < vector<double> >& features, const vector<int>& sizes)
{
    for (unsigned int i = 0; i < features.size(); i++)
        dump_feature(features[i], sizes[i]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Parse command line
//////////////////////////////////////////////////////////////////////////////////////////////////

char** parse_model_feature(char **argv, int& nf, ModelFeatures& m_features)
{
    nf = atoi(argv[0]);
    argv++;
    
    m_features.size = nf;
    m_features.fields.resize(nf);
    m_features.fns.resize(nf);
    m_features.sizes.resize(nf);
    
    const int items_per_field = 3;
    for (int i = 0; i < nf; i++)
    {    
        m_features.fields[i] = argv[i*items_per_field];
        m_features.sizes[i] = atoi(argv[i*items_per_field+1]);
        m_features.fns[i] = argv[i*items_per_field+2];
    }
    return (argv + items_per_field*nf);
}

char** parse_output_features(char **argv, int& no, OutputFeatures& obins)
{    
    no = atoi(argv[0]);
    argv++;

    obins.size = no;
    obins.fields.resize(no);
    obins.from.resize(no);
    obins.to.resize(no);
    
    const int items_per_field = 3;
    for (int i = 0; i < no; i++)
    {    
        obins.fields[i] = argv[i*items_per_field];
        obins.from[i] = atoi(argv[i*items_per_field+1]);
        obins.to[i] = atoi(argv[i*items_per_field+2]);
    }
    return (argv + items_per_field*no);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Fend list intersection
//////////////////////////////////////////////////////////////////////////////////////////////////

struct comp_fends : public binary_function<Fend, Fend, bool> {
	bool operator()(Fend x, Fend y) { return x.index < y.index; }
    };

void intersect_fends(vector<Fend>& fends_x, vector<Fend>& fends_y,
                     vector<Fend>& fends_x_unique, vector<Fend>& fends_y_unique, vector<Fend>& fends_common)
{
    sort(fends_x.begin(), fends_x.end(), comp_fends());
    sort(fends_y.begin(), fends_y.end(), comp_fends());

    // X && Y
    set_intersection(fends_x.begin(), fends_x.end(), 
                     fends_y.begin(), fends_y.end(), back_inserter(fends_common), comp_fends());

    // X \ Y
    set_difference(fends_x.begin(), fends_x.end(), 
                   fends_y.begin(), fends_y.end(), back_inserter(fends_x_unique), comp_fends());
    
    // Y \ X
    set_difference(fends_y.begin(), fends_y.end(), 
                   fends_x.begin(), fends_x.end(), back_inserter(fends_y_unique), comp_fends());
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////////////////////////////////////////////

void print_user_arguments(ModelFeatures& m_features, OutputFeatures& o_features_x, OutputFeatures& o_features_y,
                          const int from_fend, const int to_fend)
{
    if (m_features.size == 0)
        cout << "No model functions" << endl;
    for (int i=0; i<m_features.size; i++)
        cerr << "M" << i << ": " << m_features.fields[i] << ", " <<
            m_features.fns[i] << ", dim=" << m_features.sizes[i] << endl;
    for (int i=0; i<o_features_x.size; i++)
        cerr << "Output bin X" << i << ": " << o_features_x.fields[i] << ", " <<
            o_features_x.from[i] << "-" << o_features_x.to[i] << endl;
    for (int i=0; i<o_features_y.size; i++)
        cerr << "Output bin Y" << i << ": " << o_features_y.fields[i] << ", " <<
            o_features_y.from[i] << "-" << o_features_y.to[i] << endl;
    cerr << "Y axis limited: from fend " << from_fend << ", to fend " << to_fend << endl;
}    

void init_keys(const OutputFeatures& features, int index, string key, vector<string>& keys)
{
    if (index == features.size) {
        keys.push_back(key);
        return;
    }
    for (int i = features.from[index]; i <= features.to[index]; i++)
    {
        string next_key = (index == 0) ? NumberToString(i) : key + string("\t") + NumberToString(i);
        init_keys(features, index+1, next_key, keys);
    }
}

void usage(const char* name)
{
    fprintf(stderr, "usage: %s <input file> <output file> <from fend> <to fend> <filter> <cis threshold> <prior> <model features> <output features x> <output features y>\n", name);
    cerr << "Each model feature: <field size fn>" << endl;
    cerr << "Each output feature: <field from_level to_level>" << endl;
    fprintf(stderr, "For example: %s fends.table out trans 100000 0.1 1 frag_len_bin frag_len.f 50 1 gc_bin 1 50 1 gc_bin 1 50\n", name);
    
    exit(1);
}

double fends_expected(const vector<Fend>& fends_x, const vector<Fend>& fends_y,
                      Filter filter, const double cis_threshold,
                      const vector < vector<double> >& model_matrices, const vector<int>& msizes)
{
    double total = 0;
    vector<Fend>::const_iterator jt_x, jt_y;
    
    for(jt_x = fends_x.begin(); jt_x != fends_x.end(); jt_x++)
    for(jt_y = fends_y.begin(); jt_y != fends_y.end(); jt_y++)
    {
        const Fend& fend_x = *jt_x;
        const Fend& fend_y = *jt_y;
        //cout << fend_x.chrom_code << " " << fend_y.chrom_code << endl;
        switch ( filter ) {
        case f_trans:
            if (fend_x.chrom_code == fend_y.chrom_code)
                continue;
            break;
        case f_close_cis:
            if (fend_x.chrom_code != fend_y.chrom_code || (abs(fend_x.coord - fend_y.coord) > cis_threshold))
                continue;
            break;
        case f_far_cis:
            if (fend_x.chrom_code != fend_y.chrom_code || (abs(fend_x.coord - fend_y.coord) < cis_threshold))
                continue;
            break;
        case f_both:
            break;
            }

        // skip if self-interaction
        if (fend_x.index == fend_y.index)
            continue;
        
        double value = 1;
        for(unsigned int mi = 0; mi < model_matrices.size(); mi++)
            value *= model_matrices[mi][MATRIX_INDEX(fend_x.levels[mi], fend_y.levels[mi], msizes[mi])];
        
        total += value;
    }
    return total;
}

int main(int argc, char **argv)
{
    if (argc == 1)
        usage(argv[0]);

    string ifn(argv[1]);
    string ofn(argv[2]);
    int from_fend = atoi(argv[3]);
    int to_fend = atoi(argv[4]);
    string filter_str(argv[5]);
    double cis_threshold = atof(argv[6]);
    double prior = atof(argv[7]);
    argv += 8;
    cerr << "Prior: " << prior << endl;

    Filter filter;
    if (filter_str == "both")
        filter = f_both;
    else if (filter_str == "trans")
        filter = f_trans;
    else if (filter_str == "close_cis")
        filter = f_close_cis;
    else if (filter_str == "far_cis")
        filter = f_far_cis;
    else {
        cerr << "Unknown filter " << filter_str << endl;
        exit(1);
    }

    if ( ((filter == f_far_cis) || (filter == f_close_cis)) && cis_threshold <= 0) {
        cerr << "Cis threshold must be positive when filter is set to close or far cis" << endl;
        exit(1);
    }
    cerr << "Filter: " << filter_str << ", cis threshold: " << cis_threshold << endl;
    
    // parse model features
    int n_m_features;
    ModelFeatures m_features;
    argv = parse_model_feature(argv, n_m_features, m_features);
    
    // parse output bins
    int n_o_features_x, n_o_features_y;
    OutputFeatures o_features_x, o_features_y;
    argv = parse_output_features(argv, n_o_features_x, o_features_x);
    argv = parse_output_features(argv, n_o_features_y, o_features_y);

    int dim_x = 1;
    for (int i=0; i<o_features_x.size; i++)
        dim_x *= o_features_x.to[i] - o_features_x.from[i] + 1;
    int dim_y = 1;
    for (int i=0; i<o_features_y.size; i++)
        dim_y *= o_features_y.to[i] - o_features_y.from[i] + 1;
    
    // print arguments back to user
    print_user_arguments(m_features, o_features_x, o_features_y, from_fend, to_fend);

    // container of output bins
    map<string, Bin> obins_x;
    map<string, Bin> obins_y;

    // list of model features
    vector < vector<double> > model_matrices(n_m_features);

    // read fends file
    cerr << "Reading fends into memory ..." << endl;
    read_fends_table(ifn, m_features, o_features_x, o_features_y, from_fend, to_fend, obins_x, obins_y);

    // read model feature files
    read_feature_tables(m_features.fns, m_features.sizes, model_matrices);

    // dump_features(model_matrices, m_features.sizes);
    // dump_bins(obins_x, "X");
    // dump_bins(obins_y, "Y");
    
    // open output file
    ofstream out(ofn.c_str());
    massert(out.is_open(), "could not open file %s", ofn.c_str());

    out << join_strings(o_features_x.fields.begin(), o_features_x.fields.end(), "\t", NumberToString(1)) << "\t";
    out << join_strings(o_features_y.fields.begin(), o_features_y.fields.end(), "\t", NumberToString(2)) << "\t";
    out << "value\n";

    int total_bins = dim_x * dim_y;
    cerr << "Traversing all fend pairs, number of bins: " << total_bins << " bins" << endl;
    int count = 0;

    vector<string> keys_x;
    vector<string> keys_y;
    init_keys(o_features_x, 0, string(""), keys_x);
    init_keys(o_features_y, 0, string(""), keys_y);

    vector<int>& msizes = m_features.sizes;
    
    for (unsigned int ix=0; ix<keys_x.size(); ix++)
    for (unsigned int iy=0; iy<keys_y.size(); iy++)
    {
        map<string, Bin>::iterator it_x = obins_x.find(keys_x[ix]);
        map<string, Bin>::iterator it_y = obins_y.find(keys_y[iy]);

        if (it_x == obins_x.end() || it_y == obins_y.end())
        {
            out << keys_x[ix] << "\t" << keys_y[iy] << "\t" << 0 << endl;
            continue;
        }
        
        vector<Fend>& fends_x = (*it_x).second.fends;
        vector<Fend>& fends_y = (*it_y).second.fends;

        vector<Fend> fends_x_unique;
        vector<Fend> fends_y_unique;
        vector<Fend> fends_common;
        intersect_fends(fends_x, fends_y, fends_x_unique, fends_y_unique, fends_common);
        
//         double total = fends_expected(fends_x_unique, fends_y_unique, filter, cis_threshold, model_matrices, msizes) +
//                        fends_expected(fends_common, fends_common, filter, cis_threshold, model_matrices, msizes)/2;
//         double total = fends_expected(fends_x, fends_y, filter, cis_threshold, model_matrices, msizes) -
//                        fends_expected(fends_common, fends_common, filter, cis_threshold, model_matrices, msizes)/2;
        double total = fends_expected(fends_x, fends_y, filter, cis_threshold, model_matrices, msizes);
        if (keys_x[ix] == keys_y[iy])
            total = total / 2;
        
        total *= prior;

        if (total > 0)
            out << keys_x[ix] << "\t" << keys_y[iy] << "\t" << total << endl;

        // progress report
        if (++count % 100 == 0)
            fprintf(stderr, "%d/%d (%.1f%%)\n", count, total_bins, 100.0 * count / total_bins);
    }
    fprintf(stderr, "%d/%d (%.1f%%)\n", count, total_bins, 100.0 * count / total_bins);
    
    out.close();
    cerr << "done" << endl;
    
    return 0;
}
