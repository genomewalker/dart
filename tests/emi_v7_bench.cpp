// Real-data bench: parse a 16-col mmseqs tsv -> v7 emi. Measures compression
// (vs gzipped tsv, vs the current v6 emi) and write/read speed (no regression?).
#include "dart/emi_v7_io.hpp"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <fstream>
#include <sstream>
#include <unordered_map>

using namespace dart::emiv7;
using clk = std::chrono::steady_clock;
static double ms(clk::time_point a, clk::time_point b){ return std::chrono::duration<double,std::milli>(b-a).count(); }

int main(int argc, char** argv) {
    if (argc < 3) { printf("usage: bench <in.tsv> <out.emi>\n"); return 2; }
    const char* tsv = argv[1]; const char* out = argv[2];

    std::vector<Row> rows;
    std::vector<std::string> read_names, ref_names;
    std::unordered_map<std::string,uint32_t> rmap, fmap;
    auto did = [](std::unordered_map<std::string,uint32_t>& m, std::vector<std::string>& names, const std::string& s)->uint32_t{
        auto it=m.find(s); if(it!=m.end()) return it->second; uint32_t id=names.size(); m.emplace(s,id); names.push_back(s); return id; };

    auto t0 = clk::now();
    std::ifstream in(tsv);
    if (!in) { printf("cannot open %s\n", tsv); return 1; }
    std::string line; size_t bytes=0;
    while (std::getline(in, line)) {
        bytes += line.size()+1;
        // query target fident alnlen mismatch gapopen qstart qend tstart tend evalue bits qlen tlen qaln taln
        std::vector<std::string> c; c.reserve(16); size_t p=0;
        for (int f=0; f<16; ++f){ size_t q=line.find('\t',p); std::string tok=line.substr(p, q==std::string::npos?std::string::npos:q-p); c.push_back(tok); if(q==std::string::npos){ break;} p=q+1; }
        if (c.size()<16) continue;
        Row r;
        r.read_idx = did(rmap, read_names, c[0]);
        r.ref_idx  = did(fmap, ref_names, c[1]);
        r.bit_score = std::strtof(c[11].c_str(),nullptr);
        r.damage_score = 0.0f;
        double ev = std::strtod(c[10].c_str(),nullptr);
        r.evalue_log10 = (ev>0)? (float)std::log10(ev) : -300.0f;
        r.tstart=(uint16_t)std::strtoul(c[8].c_str(),nullptr,10);
        r.tend  =(uint16_t)std::strtoul(c[9].c_str(),nullptr,10);
        r.qlen  =(uint16_t)std::strtoul(c[12].c_str(),nullptr,10);
        r.tlen  =(uint16_t)std::strtoul(c[13].c_str(),nullptr,10);
        r.qaln = c[14]; r.taln = c[15];
        rows.push_back(std::move(r));
    }
    auto t1 = clk::now();
    std::vector<EmiReadAnnot> annots(read_names.size());   // damage payload: zeroed here (no .agd in bench)
    printf("parsed rows=%zu reads=%zu refs=%zu tsv=%.1fMB parse=%.0fms\n",
           rows.size(), read_names.size(), ref_names.size(), bytes/1e6, ms(t0,t1));

    auto t2 = clk::now();
    Writer::write(out, rows, annots, read_names, ref_names);
    auto t3 = clk::now();
    std::ifstream es(out, std::ios::ate|std::ios::binary); size_t esz=es.tellg();
    printf("v7 write=%.0fms  emi=%.2fMB  (%.1f%% of tsv)\n", ms(t2,t3), esz/1e6, 100.0*esz/bytes);

    // read: full hot scan (EM path) + random access sample
    auto t4 = clk::now();
    Reader rd(out);
    uint64_t nseen=0; double sbit=0;
    for (uint32_t b=0;b<rd.num_blocks();++b){ auto hb=rd.hot_block(b); nseen+=hb.n; for(uint32_t j=0;j<hb.n;++j) sbit+=hb.bit[j]; }
    auto t5 = clk::now();
    // random-access 1000 reads
    uint64_t racc=0; for(uint32_t r=0;r<rd.num_reads() && r<1000;++r){ auto g=rd.rows_for_read(r); racc+=g.size(); }
    auto t6 = clk::now();
    printf("v7 open+hotscan=%.0fms (%llu rows, chk=%.0f)  randacc1000=%.0fms (%llu rows)\n",
           ms(t4,t5),(unsigned long long)nseen,sbit, ms(t5,t6),(unsigned long long)racc);
    printf("v7 read_name[0]=%s ref_name[0]=%s\n", rd.read_name(0).c_str(), rd.ref_name(0).c_str());
    return 0;
}
