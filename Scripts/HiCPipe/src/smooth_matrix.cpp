#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////
// Data structures
//////////////////////////////////////////////////////////////////////////////////////////////////

struct Cbin
{
    string chr;
    double from, to, center;
};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Utility functions
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
// I/O functions
//////////////////////////////////////////////////////////////////////////////////////////////////

void read_cbins_table(string fn, vector< Cbin > &cbins)
{
    cerr << "reading cbins file " << fn << endl;
    
    ifstream in(fn.c_str());
    assert(in.is_open());
	vector<string> fields;
    char delim = '\t';
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
        
        int index = atoi(fields[0].c_str());
        
        Cbin cbin;
        cbin.chr = fields[1];
        cbin.from = atof(fields[2].c_str());
        cbin.to = atof(fields[3].c_str());
        cbin.center = (cbin.from+cbin.to)/2;

        // add as last element, make sure index is correct
        cbins.push_back(cbin);
		assert(cbins.size() == (unsigned int)index);
	}
    
    in.close();
}

void read_matrix_table(string fn, vector<double>& matrix, int dim,
                       string field_bin_title, string field_value_title)
{
    cerr << "reading matrix file " << fn << endl; 
    ifstream in(fn.c_str());
    assert(in.is_open());
	vector<string> fields;
    char delim = '\t';
    int matrix_size = matrix.size();

    split_line(in, fields, delim);
    int value_ind = get_field_index(field_value_title, fields);
    int bin_ind1 = get_field_index(field_bin_title + "1", fields);
    int bin_ind2 = get_field_index(field_bin_title + "2", fields);

	while(in) {
		split_line(in, fields, delim);
		if(fields.size() == 0)
			return;
        
        int i = atoi(fields[bin_ind1].c_str()) - 1;
        int j = atoi(fields[bin_ind2].c_str()) - 1;
        double value = atof(fields[value_ind].c_str());

        int index1 = i + j*dim;
        int index2 = j + i*dim;
        
        assert( (index1 >= 0) && (index1 < matrix_size) &&
                (index2 >= 0) && (index2 < matrix_size) );

        matrix[index1] = value;
        matrix[index2] = value;
	}
    
    in.close();
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////////////////////////////////////////////

void usage(const char* name)
{
    fprintf(stderr,
            "usage: %s <input cbins> <input matrix> <bin field> <value field> <step> <smooth width> <normalize T|F> <output file>\n",
            name);
    exit(1);
}

int main(int argc, char **argv)
{
    if (argc == 1)
        usage(argv[0]);

    assert(argc == 9);
    string cbins_fn(argv[1]);
    string matrix_fn(argv[2]);
    string field_bin_title(argv[3]);
    string field_value_title(argv[4]);
    int step = atoi(argv[5]);
    int width = atoi(argv[6]);
    string normalize_str(argv[7]);
    string ofn(argv[8]);

    assert(normalize_str == "T" || normalize_str == "F");
    bool normalize = (normalize_str == "T");
    
    // read cbins
    vector<Cbin> cbins;
    read_cbins_table(cbins_fn, cbins);
    int dim = cbins.size();
        
    // read frag_len function
    int matrix_size = dim*dim;
    vector<double> matrix(matrix_size);
    read_matrix_table(matrix_fn, matrix, dim, field_bin_title, field_value_title);

    // open output file
    cerr << "Creating output " << ofn << endl;
    ofstream out(ofn.c_str());
    assert(out.is_open());
    out << "cbin1\tchr1\tfrom1\tto1\tcbin2\tchr2\tfrom2\tto2\t" << field_value_title << endl;

    int count = 0;
    cerr << "smoothing matrix of dimension " << dim << endl;
    double max_total_weight = 0;
    
    for (int i1 = 0; i1 < dim; i1++)
    for (int i2 = 0; i2 < dim; i2++)
    {    
        // progress report
        if (++count % 100000 == 0)
            fprintf(stderr, "%d/%d (%.1f%%)\n", count, matrix_size, 100.0 * count / matrix_size);
        
        string chr1 = cbins[i1].chr;
        double from1 = cbins[i1].from;
        double to1 = cbins[i1].to;
        double coord1 = cbins[i1].center;
        string chr2 = cbins[i2].chr;
        double from2 = cbins[i2].from;
        double to2 = cbins[i2].to;
        double coord2 = cbins[i2].center;
        
        double total_value = 0;
        double total_weight = 0;
        int contrib_count = 0;
        
        // cerr << "********* ("<< i1 << "," << i2 << ")" << endl;
        for (int j1 = (i1-width); j1 <= (i1+width); j1++)
        for (int j2 = (i2-width); j2 <= (i2+width); j2++)
        {
            if ( (j1 < 0) || (j1 >= dim) || (j2 < 0) || (j2 >= dim) ||
                 (cbins[j1].chr != chr1) || (cbins[j2].chr != chr2) )
                continue;

            double dist1 = (coord1 - cbins[j1].center) / step;
            double dist2 = (coord2 - cbins[j2].center) / step;
            
            // Use L1 metric
            double dist = fabs(dist1) + fabs(dist2);
            
            double weight = 1  - dist / (width + 1);
            if (weight <= 0)
                continue;
            //cerr << "dist: " << dist << ", weight:" << weight << endl;
            
            int index = j1 + j2*dim;
            assert( (index >= 0) && (index < matrix_size) );
            double value = matrix[index];
            //cerr << "("<< j1 << "," << j2 << "):" << ", v=" << value << ", w=" << weight << endl;
            total_weight += weight;
            total_value += value*weight;
            contrib_count++;
        }
        double mean_value = 0;
        if (contrib_count > 0)
        {
            assert(total_weight != 0);
            mean_value = total_value;
            if (normalize)
                mean_value /= total_weight;
        }
        // cerr << "mean: " << total_value << "/" << total_weight << "=" << mean_value << endl;
        // cerr << "contrib: " << contrib_count << endl;
        // cerr << "total weights: " << total_weight << endl;
        if (max_total_weight < total_weight)
            max_total_weight = total_weight;
        // output 1-based indices
        out << i1+1 << "\t" << chr1 << "\t" << from1 << "\t" << to1 << "\t"
            << i2+1 << "\t" << chr2 << "\t" << from2 << "\t" << to2 << "\t"
            << mean_value << endl;
    }
    fprintf(stderr, "%d/%d (%.1f%%)\n", count, matrix_size, 100.0 * count / matrix_size);
    out.close();
    // cerr << "done" << endl;
    cerr << "done. max total weights: " << max_total_weight << endl;
    return 0;
}
