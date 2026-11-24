// Bounded-memory gzip scanner. Python is the reference implementation.
#include <zlib.h>
#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using U = std::uint64_t;

int main(int argc, char** argv) {
    try {
        if (argc < 3) throw std::runtime_error("Expected checkpoint file and raw files");
        std::ifstream grid(argv[1]);
        std::string prime_line;
        std::getline(grid, prime_line);
        std::istringstream prime_input(prime_line);
        std::vector<U> primes;
        U v;
        while (prime_input >> v) {
            if (v != 3 && v != 5 && v != 7 && v != 11 && v != 13 && v != 17 && v != 19)
                throw std::runtime_error("Unsupported prime");
            primes.push_back(v);
        }
        if (primes.empty()) throw std::runtime_error("No primes selected");
        std::vector<U> xs;
        while (grid >> v) xs.push_back(v);
        if (xs.empty()) throw std::runtime_error("Empty checkpoint grid");
        using Histogram = std::map<std::vector<U>, U>;
        std::vector buckets(primes.size(), std::vector<Histogram>(xs.size()));
        std::vector<U> values, parts;
        const std::regex name("(?:^|/)cl([0-9]+)mod([0-9]+)\\.([0-9]+)\\.gz$");
        for (int file=2; file<argc; ++file) {
            std::string path(argv[file]);
            std::smatch match;
            if (!std::regex_search(path, match, name)) throw std::runtime_error("Invalid raw filename");
            U r=std::stoull(match[1]), m=std::stoull(match[2]);
            U k=std::stoull(match[3]);
            if (k>=4096 || !((m==8 && (r==3 || r==7)) || (m==16 && (r==4 || r==8))))
                throw std::runtime_error("Invalid raw filename");
            const U width=1ULL<<28, stop=(k+1)*width;
            U D=k*width+r, previous=0, count=0;
            size_t index=0;
            gzFile in=gzopen(path.c_str(), "rb");
            if (!in) throw std::runtime_error("Cannot open "+path);
            gzbuffer(in, 1<<20);
            char line[16384];
            while (gzgets(in, line, sizeof(line))) {
                const auto size = std::strlen(line);
                if (size == sizeof(line)-1 && line[size-1] != '\n')
                    throw std::runtime_error("Overlong raw record");
                values.clear();
                const char *p = line, *end = line + size;
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
                if (D<=previous || D>=stop)
                    throw std::runtime_error("Unsorted or out-of-range discriminant in "+path);
                previous=D;
                while (index<xs.size() && xs[index]<D) ++index;
                if (index==xs.size()) break;
                U product = 1;
                for (size_t j = 2; j < values.size(); ++j) {
                    const U n = values[j];
                    if (!n || product > std::numeric_limits<U>::max() / n)
                        throw std::runtime_error("Invalid invariant factor");
                    product *= n;
                }
                if (product != values[1])
                    throw std::runtime_error("Class number does not match invariants");
                for (size_t t = 0; t < primes.size(); ++t) {
                    const U prime = primes[t];
                    parts.clear();
                    if (product % prime == 0) {
                        for (size_t j = 2; j < values.size(); ++j) {
                            U n = values[j], q = 1;
                            while (n % prime == 0) {
                                q *= prime;
                                n /= prime;
                            }
                            if (q > 1) parts.push_back(q);
                        }
                        std::sort(parts.rbegin(), parts.rend());
                    }
                    ++buckets[t][index][parts];
                }
                ++count;
            }
            int err=Z_OK;
            const char* msg=gzerror(in,&err);
            std::string error=msg?msg:"gzip error";
            int closed=gzclose(in);
            if ((err!=Z_OK && err!=Z_STREAM_END) || closed!=Z_OK)
                throw std::runtime_error(error);
            std::cerr<<path<<": "<<count<<" records\n";
        }
        for (size_t t = 0; t < primes.size(); ++t) {
            for (size_t i = 0; i < xs.size(); ++i) {
                for (const auto& [parts, count] : buckets[t][i]) {
                    std::string key;
                    for (U q : parts) {
                        if (!key.empty()) key += "x";
                        key += "C" + std::to_string(q);
                    }
                    if (key.empty()) key = "trivial";
                    std::cout << primes[t] << ' ' << i << ' ' << key << ' ' << count << '\n';
                }
            }
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n'; return 1;
    }
}
