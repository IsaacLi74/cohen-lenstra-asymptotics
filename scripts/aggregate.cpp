// Bounded-memory gzip scanner. Python is the reference implementation.
#include <zlib.h>
#include <algorithm>
#include <charconv>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

using U = std::uint64_t;

int main(int argc, char** argv) {
    try {
        if (argc < 3) throw std::runtime_error("Expected checkpoint file and raw files");
        std::ifstream grid(argv[1]);
        std::vector<U> xs;
        U v;
        while (grid >> v) xs.push_back(v);
        if (xs.empty()) throw std::runtime_error("Empty checkpoint grid");
        std::vector<std::map<std::string,U>> buckets(xs.size());
        const std::regex name("(?:^|/)cl([0-9]+)mod([0-9]+)\\.0\\.gz$");
        for (int file=2; file<argc; ++file) {
            std::string path(argv[file]);
            std::smatch match;
            if (!std::regex_search(path, match, name)) throw std::runtime_error("Invalid raw filename");
            U r=std::stoull(match[1]), m=std::stoull(match[2]);
            if (!((m==8 && (r==3 || r==7)) || (m==16 && (r==4 || r==8))))
                throw std::runtime_error("Invalid raw filename");
            U D=r, previous=0, count=0;
            size_t index=0;
            gzFile in=gzopen(path.c_str(), "rb");
            if (!in) throw std::runtime_error("Cannot open "+path);
            gzbuffer(in, 1<<20);
            char line[16384];
            while (gzgets(in, line, sizeof(line))) {
                std::string text(line);
                if (text.size()==sizeof(line)-1 && text.back()!='\n')
                    throw std::runtime_error("Overlong raw record");
                std::vector<U> values;
                const char *p=text.data(), *end=p+text.size();
                while (p<end) {
                    while (p<end && (*p==' ' || *p=='\t' || *p=='\r' || *p=='\n')) ++p;
                    if (p==end) break;
                    auto parsed=std::from_chars(p,end,v);
                    if (parsed.ec!=std::errc()) throw std::runtime_error("Invalid integer in "+path);
                    values.push_back(v); p=parsed.ptr;
                }
                if (values.size()<3) throw std::runtime_error("Malformed record in "+path);
                if (values[0]>(std::numeric_limits<U>::max()-D)/m)
                    throw std::runtime_error("Discriminant overflow");
                D+=m*values[0];
                if (D<=previous || D>=(1ULL<<28))
                    throw std::runtime_error("Unsorted or out-of-range discriminant in "+path);
                previous=D;
                while (index<xs.size() && xs[index]<D) ++index;
                if (index==xs.size()) break;
                U product=1;
                std::vector<U> parts;
                for (size_t j=2;j<values.size();++j) {
                    U n=values[j], q=1;
                    if (!n || product>std::numeric_limits<U>::max()/n)
                        throw std::runtime_error("Invalid invariant factor");
                    product*=n;
                    while (n%3==0) { q*=3; n/=3; }
                    if (q>1) parts.push_back(q);
                }
                if (product!=values[1] || (values[1]%3==0)!=!parts.empty())
                    throw std::runtime_error("Class number does not match invariants");
                std::sort(parts.rbegin(),parts.rend());
                std::string key;
                for (U q:parts) { if (!key.empty()) key+="x"; key+="C"+std::to_string(q); }
                if (key.empty()) key="trivial";
                ++buckets[index][key]; ++count;
            }
            int err=Z_OK;
            const char* msg=gzerror(in,&err);
            std::string error=msg?msg:"gzip error";
            int closed=gzclose(in);
            if ((err!=Z_OK && err!=Z_STREAM_END) || closed!=Z_OK)
                throw std::runtime_error(error);
            std::cerr<<path<<": "<<count<<" records\n";
        }
        for (size_t i=0;i<xs.size();++i)
            for (const auto& entry:buckets[i])
                std::cout<<i<<' '<<entry.first<<' '<<entry.second<<'\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n'; return 1;
    }
}
