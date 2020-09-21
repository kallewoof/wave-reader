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
    uint16_t format;            // 1 = PCM
    uint16_t channels;          // 1 = mono, 2 = stereo
    uint32_t sample_rate;       // samples per second, e.g. 44100
    uint32_t bytes_per_sec;     // sample rate * block align
    uint16_t block_align;       // channels * bits/sample / 8
    uint16_t bits_per_sample;   // 8 or 16
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
            // apparently 1-byte samples are unsigned; we undo that
            uint8_t i8;
            memcpy(&i8, &pb->riff->data->data[pb->pos++], 1);
            pb->chn[channel] = (int16_t)i8 - 127;
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
    struct analyzer_t* anl = malloc(sizeof(struct analyzer_t));
    anl->pb = pb;
    anl->minrange = minrange;
    anl->maxrange = maxrange;
    anl->topcount = anl->next = anl->nextnext = 0;
    anl->steps = steps;
    anl->counts = malloc(steps * sizeof(uint32_t));
    return anl;
}

bool analyzer_analyze(struct analyzer_t* anl, size_t samples) {
    memset(anl->counts, 0, anl->steps * sizeof(uint32_t));
    int16_t div = (anl->maxrange - anl->minrange + anl->steps) / anl->steps;
    div += !div;
    anl->topcount = anl->next = anl->nextnext = 0;
    for (size_t i = 0; i < samples; ++i) {
        if (!playback_iterate(anl->pb)) return false;
        int32_t c = anl->pb->chn[0];
        c -= anl->minrange;
        c /= div;
        if (c < 0 || c >= anl->steps) printf("(%d - %d) / %d = %d (outside 0..%u)\n", anl->pb->chn[0], anl->minrange, div, c, anl->steps);
        assert(c > -1 && c < anl->steps);
        if (++anl->counts[c] > anl->counts[anl->topcount]) {
            if (anl->next != c) anl->nextnext = anl->next;
            anl->next = anl->topcount;
            anl->topcount = c;
        }
    }
    return true;
}

const char* sample_timestr(uint32_t sample, uint32_t rate) {
    static char* rv = NULL;
    if (!rv) rv = malloc(6);
    uint32_t hs = sample * 100 / rate;
    sprintf(rv, "%02u.%02u", hs / 100, hs % 100);
    return rv;
}

int main(const int argc, const char* argv[]) {
    /*
     * 32-64 bits of timestamp
     * 16-32 bits of temp (0000..9999)
     * 7 (days)
     * = 336 .. 672 bits of data
     * transmitted over 29.6 seconds with 44100 sample rate
     * total 1305771 samples; given pauses in between each burst of data, assume 80% was transmission and 20% was gap
     * giving 23.7 seconds of transmission to transmit 7680-30720 bits of data: 0.0007714844 up to 0.003085937 seconds per bit
     * i.e. 40-200 samples per bit
     *
     * each transmission burst has ~52000 samples over ~1.2 seconds.
     * the data is split into 15 bursts. if we assume no repetition of data,
     * and assume that the data is simply cut into 15 pieces arbitrarily, we get something like this:
     *        |       |       |       |       |       |       |       |       |       |       |       |       |       |       
     * [timestamp][temp][timestamp][temp][timestamp][temp][timestamp][temp][timestamp][temp][timestamp][temp][timestamp][temp]
     * which could be the case e.g. if the data is serialized into a data stream that is abstracted away from the actual
     * transmission details.
     * if however each burst is a day, there is some repetition and/or checksumming as we are receiving 2x + 1 days' worth of
     * data
     * it is also possible that the data is compressed before being transmitted. that sounds overkill tho.
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

    playback_rewind(pb);

    // determine "making sound" period and "not making sound" period; there should be 15 of them
    // (and bonus also detect the bleeps)

    int making_sound = 0;
    int silent = 0;
    uint32_t sample = 0;

    /*
     * sample       time    delta   event
     * 64112        01.454          initial beep begin
     * 69689        01.580  00.126  end
     * 104486       02.369  00.915  large click begin (delta very close to 1 s from initial beep begin)
     * 104688       02.374  00.005  end
     * 108671       02.464  01.010  after-click begin (like a button coming back up after being pressed; very very close to 1s from in)
     * 109400       02.481  00.017  end
     * 118125       02.679          first transmission begin
     * 170316       03.862  01.183  end (52191 samples)
     * 174713       03.962  00.077  after?-click begin (delta is for time after first trans end)
     * 175623       03.982  00.020  end
     * 179369       04.067  00.182  after??-click begin (delta is fo rtime after first trans end)
     * 180007       04.082  00.015  end
     * 184266       04.178          second transmission begin
     * 236186       05.356  01.178  end (51920 samples)
     * 
     */
    #define SWITCH_THRESH   10
    #define SOUND_THRESH    150
    #define SILENCE_THRESH  100
    while (playback_iterate(pb)) {
        ++sample;
        int16_t c = pb->chn[0];
        if (abs(c) > SOUND_THRESH) {
            ++making_sound;
            if (making_sound == SWITCH_THRESH) {
                printf("%-13u %s: making sound\n", sample, sample_timestr(sample, riff->fmt->sample_rate));
            }
            if (making_sound >= SWITCH_THRESH && making_sound < SWITCH_THRESH + 10) {
                printf("%-13u %s %u %d %d\n", sample, sample_timestr(sample, riff->fmt->sample_rate), abs(c), making_sound, silent);
            }
            if (making_sound >= SWITCH_THRESH && silent > 0) if (silent >= SWITCH_THRESH) silent = 0; else --silent;
        } else if (abs(c) < SILENCE_THRESH) {
            ++silent;
            // if (silent < SWITCH_THRESH) {
            //     printf("%s %u %d %d\n", sample_timestr(sample, riff->fmt->sample_rate), abs(c), making_sound, silent);
            // }
            if (silent == SWITCH_THRESH) {
                printf("%-13u %s: silent\n", sample, sample_timestr(sample, riff->fmt->sample_rate));
            }
            if (silent >= SWITCH_THRESH && making_sound > 0) if (making_sound >= SWITCH_THRESH) making_sound = 0; else --making_sound;
        }
    }

    // // we define quality of parameters as the set where we saw the highest top count divided by total counts,
    // // i.e. a perfect score for 100 samples is where top sample count = 100
    // // and a 50% score is where the top sample count = 50, so we can calculate the quality
    // // of a single analyze as topcount / count,
    // // and the aggregate score as the sum of all topcounts divided by the sum of counts
    // // our goal is to find a signal with ON and OFF (true/false, 1/0) elements, which means we also want
    // // to raise quality for when there are two distinct signals, where the distinction is defined
    // // by the steps implementation (i.e. a 58 is considered distinctly different from a 59)
    // // finally, the end quality is defined as the top TWO hits across the entire run, and their counts over the entire counts
    // double top_score = 0;
    // #define STEPS 1000
    // struct analyzer_t* anl = analyzer_init(pb, min, max, STEPS);
    // uint64_t tally_counts[STEPS];
    // // for (size_t window = 30; window < 500; ++window) {
    // //     double local_score = 0;
    // //     size_t local_slide = 0;
    // //     uint32_t local_top[4] = {0,0,0,0};
    // //     uint64_t local_tally_sum = 0;
    // //     for (size_t slide = 0; slide < window; ++slide) {
    // //         printf("w = %5zu, slide = %5zu\r", window, slide);
    // //         fflush(stdout);
    //         playback_rewind(pb);
    // //         // "slide"
    // //         for (size_t i = 0; i < slide; ++i) playback_iterate(pb);
    // //         memset(tally_counts, 0, sizeof(uint64_t) * STEPS);
    // //         uint32_t top[4] = {0, 0, 0, 0};
    // //         uint32_t samples = 0;
    // //         while (analyzer_analyze(anl, window)) {
    // //             samples += window;
    // //             int16_t tc = anl->topcount;
    // //             tally_counts[tc] += anl->counts[tc];
    // //             // printf("%10d [%3u]\t%10d [%3u]\t%10d [%3u]\n", anl->topcount, anl->counts[anl->topcount], anl->next, anl->counts[anl->next], anl->nextnext, anl->counts[anl->nextnext]);
    // //         }
    // //         uint32_t min = 999, max = 0;
    // //         for (size_t j = 0; j < 4; ++j) {
    // //             for (size_t i = 0; i < STEPS; ++i) {
    // //                 if (tally_counts[i] > tally_counts[top[j]] && (j == 0 || tally_counts[top[j-1]] > tally_counts[i])) {
    // //                     top[j] = i;
    // //                     if (min > i) min = i;
    // //                     if (max < i) max = i;
    // //                 }
    // //             }
    // //         }
    // //         double score = (double)(tally_counts[top[0]] + tally_counts[top[1]] + tally_counts[top[2]] + tally_counts[top[3]]) * (max - min) / samples;
    // //         if (score > top_score) {
    // //             printf("new top score %16lf: window = %5zu, slide = %5zu, top = %4u, %4u, %4u, %4u [%5llu]\n", score, window, slide, top[0], top[1], top[2], top[3], tally_counts[top[0]] + tally_counts[top[1]] + tally_counts[top[2]] + tally_counts[top[3]]);
    // //             local_score = top_score = score;
    // //         }
    // //         if (score > local_score) {
    // //             local_score = score;
    // //             local_slide = slide;
    // //             memcpy(local_top, top, sizeof(uint32_t) * 4);
    // //             local_tally_sum = tally_counts[top[0]] + tally_counts[top[1]] + tally_counts[top[2]] + tally_counts[top[3]];
    // //         }
    // //     }
    // //     if (local_score < top_score) printf("  local score %16lf: window = %5zu, slide = %5zu, top = %4u, %4u, %4u, %4u [%5llu]\n", local_score, window, local_slide, local_top[0], local_top[1], local_top[2], local_top[3], local_tally_sum);
    // // }
    // uint32_t sample = 0;
    // while (playback_iterate(pb)) {
    //     sample++;
    //     uint32_t hs = sample / 441;
    //     uint32_t s = 0;
    //     if (hs > 99) {
    //         s = hs / 100;
    //         hs = hs % 100;
    //     }
    //     printf("%02u.%02u %10d %10d\n", s, hs, pb->chn[0], pb->chn[1]);
    // }
}
