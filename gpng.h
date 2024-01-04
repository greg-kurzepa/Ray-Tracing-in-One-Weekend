#ifndef GPNG
#define GPNG

#include <vector>
#include <cstdint>
#include <string>

namespace gpng {

    class Image {
        public:
            uint8_t* image;
            std::vector<uint8_t> main_buffer;

            int width;
            int height;

            Image(int w, int h);
            ~Image();

            uint8_t& operator()(int column, int row, int colour);
            void close();
            void save(std::string filename);
            void deflate_no_compression(std::vector<uint8_t>& buffer);
    };

}

#endif // GPNG