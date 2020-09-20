#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>

struct tag_t {
    char header[4];
    uint32_t len;
};

#define set_tag_header(th, h) do { th[0] = h[0]; th[1] = h[1]; th[2] = h[2]; th[3] = h[3]; } while (0)

#define fpull(stm, ob) fread(&ob, sizeof(ob), 1, stm)

struct fmt_t {
    struct tag_t tag;
    uint16_t format; // 1 = PCM
    uint16_t channels; // 1 = mono, 2 = stereo
    uint32_t sample_rate; // samples per second, e.g. 44100
    uint32_t bytes_per_sec; // sample rate * block align
    uint16_t block_align; // channels * bits/sample / 8
    uint16_t bits_per_sample; // 8 or 16
};

struct fmt_t* fmt_read(FILE* stream) {
    struct fmt_t* fmt = malloc(sizeof(struct fmt_t));
    set_tag_header(fmt->tag.header, "fmt ");
    fpull(stream, fmt->tag.len);                    assert(fmt->tag.len == 16);
    fpull(stream, fmt->format);                     assert(fmt->format == 1);
    fpull(stream, fmt->channels);                   assert(fmt->channels > 0 && fmt->channels < 3);
    fpull(stream, fmt->sample_rate);
    fpull(stream, fmt->bytes_per_sec);
    fpull(stream, fmt->block_align);
    fpull(stream, fmt->bits_per_sample);            assert(fmt->bits_per_sample == 8 || fmt->bits_per_sample == 16);
    return fmt;
}

struct data_t {
    struct tag_t tag;
    uint8_t* data;
};

// struct data_t* data_read(FILE* stream, struct data_t* data_in) {
//     struct data_t* data = data_in ?: malloc(sizeof(struct data_t));
//     uint32_t new_len;
//     fpull(stream, new_len);
//     if (data_in) {
//         data->data = realloc(data->data, data->tag.len + new_len);

//         in_len = data->tag.len;
//     } else {
//         set_tag_header(data->tag.header, "data");
//         in_len = 0;
//         data->data = NULL;
//     }
    
// }

struct data_t* data_read(FILE* stream) {
    struct data_t* data = malloc(sizeof(struct data_t));
    set_tag_header(data->tag.header, "data");
    fpull(stream, data->tag.len);
    data->data = malloc(data->tag.len);
    size_t i = fread(data->data, 1, data->tag.len, stream);
    if (i != data->tag.len) fprintf(stderr, "unable to read %u bytes of data (read %zu)\n", data->tag.len, i);
    return data;
}

struct riff_t {
    struct fmt_t* fmt;
    struct data_t* data;
};

struct riff_t* riff_read(FILE* stream)
{
    struct riff_t* riff = malloc(sizeof(struct riff_t));
    riff->fmt = NULL;
    riff->data = NULL;
    char header[5];
    header[4] = 0;
    uint32_t u32;
    while (4 == fread(header, 1, 4, stream)) {
        printf("tag: %s\n", header);
        #define eq(tag) (header[0] == tag[0] && header[1] == tag[1] && header[2] == tag[2] && header[3] == tag[3])
        if (eq("RIFF")) {
            // move beyond header cruft
            fpull(stream, u32); // file length
            fread(header, 1, 4, stream); // "WAVE"
            assert(eq("WAVE"));
        } else if (eq("fmt ")) {
            assert(NULL == riff->fmt);
            riff->fmt = fmt_read(stream);
        } else if (eq("data")) {
            assert(NULL == riff->data);
            riff->data = data_read(stream);
        } else {
            // skip unknown
            fpull(stream, u32); // chunk length
            fseek(stream, u32, SEEK_CUR);
        }
    }
    return riff;
}

struct playback_t {
    struct riff_t* riff;
    size_t pos;
    int16_t chn[2];
};

struct playback_t* playback_init(struct riff_t* riff) {
    struct playback_t* pb = malloc(sizeof(struct playback_t));
    pb->riff = riff;
    pb->pos = 0;
    return pb;
}

bool playback_iterate(struct playback_t* pb) {
    uint8_t bytes_per_sample = pb->riff->fmt->bits_per_sample >> 3;
    if (pb->pos + bytes_per_sample * pb->riff->fmt->channels > pb->riff->data->tag.len) return false;
    for (uint8_t channel = 0; channel < pb->riff->fmt->channels; ++channel) {
        if (bytes_per_sample == 1) {
            int8_t i8;
            memcpy(&i8, &pb->riff->data->data[pb->pos++], 1);
            pb->chn[channel] = i8;
        } else {
            memcpy(&pb->chn[channel], &pb->riff->data->data[pb->pos], 2);
            pb->pos += 2;
        }
    }
    return true;
}

void playback_rewind(struct playback_t* pb) {
    pb->pos = 0;
}

struct analyzer_t {
    struct playback_t* pb;
    int16_t minrange, maxrange, topcount, next, nextnext;
    uint32_t steps;
    uint32_t* counts;
};

struct analyzer_t* analyzer_init(struct playback_t* pb, int16_t minrange, int16_t maxrange, uint32_t steps) {
    struct analyzer_t* anal = malloc(sizeof(struct analyzer_t));
    anal->pb = pb;
    anal->minrange = minrange;
    anal->maxrange = maxrange;
    anal->topcount = anal->next = anal->nextnext = 0;
    anal->steps = steps;
    anal->counts = malloc(steps * sizeof(uint32_t));
    return anal;
}

bool analyzer_analyze(struct analyzer_t* anal, size_t samples) {
    memset(anal->counts, 0, anal->steps * sizeof(uint32_t));
    int16_t div = (anal->maxrange - anal->minrange + anal->steps) / anal->steps;
    div += !div;
    anal->topcount = anal->next = anal->nextnext = 0;
    for (size_t i = 0; i < samples; ++i) {
        if (!playback_iterate(anal->pb)) return false;
        int32_t c = anal->pb->chn[0];
        c -= anal->minrange;
        c /= div;
        if (c < 0 || c >= anal->steps) printf("(%d - %d) / %d = %d (outside 0..%u)\n", anal->pb->chn[0], anal->minrange, div, c, anal->steps);
        assert(c > -1 && c < anal->steps);
        if (++anal->counts[c] > anal->counts[anal->topcount]) {
            if (anal->next != c) anal->nextnext = anal->next;
            anal->next = anal->topcount;
            anal->topcount = c;
        }
    }
    return true;
}

int main(const int argc, const char* argv[]) {
    /*
     * 32-64 bits of timestamp
     * 16-32 bits of temp (0000..9999)
     * 15 (days)
     * = 7680 .. 30720 bits of data
     * transmitted over 29.6 seconds with 44100 sample rate
     * total 1305771 samples; given pauses in between each burst of data, assume 80% was transmission and 20% was gap
     * giving 23.7 seconds of transmission to transmit 7680-30720 bits of data: 0.0007714844 up to 0.003085937 seconds per bit
     * i.e. 40-200 samples per bit
     */
    if (argc != 2) {
        fprintf(stderr, "syntax: %s <wav file>\n", argv[0]);
        exit(1);
    }
    FILE* fp = fopen(argv[1], "rb");
    if (!fp) {
        fprintf(stderr, "unable to open file for reading: %s\n", argv[1]);
        exit(1);
    }
    struct riff_t* riff = riff_read(fp);
    printf("file: %s\n", argv[1]);
    printf("format: %u\n", riff->fmt->format);
    printf("channels: %u\n", riff->fmt->channels);
    printf("sample rate: %u\n", riff->fmt->sample_rate);
    printf("bytes/sec: %u\n", riff->fmt->bytes_per_sec);
    printf("block align: %u\n", riff->fmt->block_align);
    printf("bits/sample: %u\n", riff->fmt->bits_per_sample);
    printf("total data size: %u\n", riff->data->tag.len);
    fclose(fp);
    struct playback_t* pb = playback_init(riff);
    int16_t min = 99;
    int16_t max = -99;
    while (playback_iterate(pb)) {
        int16_t c = pb->chn[0];
        if (c < min) min = c;
        if (c > max) max = c;
    }
    printf("sample range: %d .. %d\n", min, max);
    // we define quality of parameters as the set where we saw the highest top count divided by total counts,
    // i.e. a perfect score for 100 samples is where top sample count = 100
    // and a 50% score is where the top sample count = 50, so we can calculate the quality
    // of a single analyze as topcount / count,
    // and the aggregate score as the sum of all topcounts divided by the sum of counts
    // our goal is to find a signal with ON and OFF (true/false, 1/0) elements, which means we also want
    // to raise quality for when there are two distinct signals, where the distinction is defined
    // by the steps implementation (i.e. a 58 is considered distinctly different from a 59)
    // finally, the end quality is defined as the top TWO hits across the entire run, and their counts over the entire counts
    double top_score = 0;
    #define STEPS 10000
    struct analyzer_t* anal = analyzer_init(pb, min, max, STEPS);
    uint64_t tally_counts[STEPS];
    for (size_t window = 10; window < 500; ++window) {
        for (size_t slide = 0; slide < window; ++slide) {
            printf("w = %5zu, slide = %5zu\r", window, slide);
            fflush(stdout);
            playback_rewind(pb);
            // "slide"
            for (size_t i = 0; i < slide; ++i) playback_iterate(pb);
            memset(tally_counts, 0, sizeof(uint64_t) * STEPS);
            uint32_t top[4] = {0, 0, 0, 0};
            uint32_t samples = 0;
            while (analyzer_analyze(anal, window)) {
                samples += window;
                int16_t tc = anal->topcount;
                tally_counts[tc] += anal->counts[tc];
                // printf("%10d [%3u]\t%10d [%3u]\t%10d [%3u]\n", anal->topcount, anal->counts[anal->topcount], anal->next, anal->counts[anal->next], anal->nextnext, anal->counts[anal->nextnext]);
            }
            uint32_t min = 999, max = 0;
            for (size_t j = 0; j < 4; ++j) {
                for (size_t i = 0; i < STEPS; ++i) {
                    if (tally_counts[i] > tally_counts[top[j]] && (j == 0 || tally_counts[top[j-1]] > tally_counts[i])) {
                        top[j] = i;
                        if (min > i) min = i;
                        if (max < i) max = i;
                    }
                }
            }
            double score = (double)(tally_counts[top[0]] + tally_counts[top[1]] + tally_counts[top[2]] + tally_counts[top[3]]) * (max - min) / samples;
            if (score > top_score) {
                printf("new top score %16lf: window = %5zu, slide = %5zu, top = %4u, %4u, %4u, %4u [%5llu]\n", score, window, slide, top[0], top[1], top[2], top[3], tally_counts[top[0]] + tally_counts[top[1]] + tally_counts[top[2]] + tally_counts[top[3]]);
                top_score = score;
            }
        }
    }
    // if (riff->fmt->channels == 2) {
    //     while (playback_iterate(pb)) {
    //         printf("%10d %10d\n", pb->chn[0], pb->chn[1]);
    //     }
    // } else {
    //     while (playback_iterate(pb)) {
    //         printf("%10d\n", pb->chn[0]);
    //     }
    // }
}
