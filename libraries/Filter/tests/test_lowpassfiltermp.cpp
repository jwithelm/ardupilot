#include <AP_gtest.h>
//#include <../../googlemock/include/gmock/gmock.h>

#include <Filter/LowPassFilter2p.h>

const AP_HAL::HAL& hal = AP_HAL::get_HAL();

typedef LowPassFilterMpFloat LPF;


/*
InitTest:
- check that sample_freq and cutoff_freq are initialised and set correctly
*/
TEST(LowPassFilterMpTest, InitTest)
{
    static constexpr float sample_freq {1000};
    static constexpr float cutoff_freq {20};
    static constexpr uint8_t f_order {3};
    static constexpr uint8_t f_types[3] {LPF::FilterType::PTn, LPF::FilterType::Bessel, LPF::FilterType::Butterworth};
    
    static constexpr uint8_t sz = ARRAY_SIZE(f_types);

    LPF lpfs[sz];
    for (uint8_t idx = 0; idx < sz; ++idx) {
        lpfs[idx].init(f_order, f_types[idx]);

        EXPECT_FLOAT_EQ(lpfs[idx].get_sample_freq(), 0.0F);
        EXPECT_FLOAT_EQ(lpfs[idx].get_cutoff_freq(), 0.0F);

        lpfs[idx].set_cutoff_frequency(sample_freq, cutoff_freq);

        EXPECT_FLOAT_EQ(lpfs[idx].get_sample_freq(), sample_freq);
        EXPECT_FLOAT_EQ(lpfs[idx].get_cutoff_freq(), cutoff_freq);
    }
}


/*
PassthroughTest:
- if either sample_freq or cutoff_freq are set to zero, the filter should just passthrough values
*/
TEST(LowPassFilterMpTest, PassthroughTest)
{
    static constexpr float sample_freq {1000};
    static constexpr float cutoff_freq {20};
    static constexpr uint8_t f_order {3};
    static constexpr uint8_t f_types[3] {LPF::FilterType::PTn, LPF::FilterType::Bessel, LPF::FilterType::Butterworth};
    
    static constexpr uint8_t sz = ARRAY_SIZE(f_types);

    float y;
    float u;

    LPF lpfs[sz];
    for (uint8_t idx = 0; idx < sz; ++idx) {
        lpfs[idx].init(f_order, f_types[idx]);

        // sample_freq & cutoff_freq not set (passthrough)
        for (uint32_t i=0; i<100; ++i) {
            u = static_cast<float>(i);
            y = lpfs[idx].apply(u);
            EXPECT_FLOAT_EQ(u, y);
        }

        // sample_freq = cutoff_freq = 0 (passthrough)
        lpfs[idx].set_cutoff_frequency(0.0F, 0.0F);
        for (uint32_t i=0; i<100; ++i) {
            u = static_cast<float>(i);
            y = lpfs[idx].apply(u);
            EXPECT_FLOAT_EQ(u, y);
        }

        // sample_freq = 1000, cutoff_freq = 0 (passthrough)
        lpfs[idx].set_cutoff_frequency(sample_freq, 0.0F);
        for (uint32_t i=0; i<100; ++i) {
            u = static_cast<float>(i);
            y = lpfs[idx].apply(u);
            EXPECT_FLOAT_EQ(u, y);
        }

        // sample_freq = 0, cutoff_freq = 20 (passthrough)
        lpfs[idx].set_cutoff_frequency(0.0F, cutoff_freq);
        for (uint32_t i=0; i<100; ++i) {
            u = static_cast<float>(i);
            y = lpfs[idx].apply(u);
            EXPECT_FLOAT_EQ(u, y);
        }

        // sample_freq = 1000, cutoff_freq = 20 (normal) -> This is defacto not a passthrough test, but a non-passthrough (filtering) test
        lpfs[idx].set_cutoff_frequency(sample_freq, cutoff_freq);
        y = lpfs[idx].apply(0);
        for (uint32_t i=1; i<100; ++i) {
            u = static_cast<float>(i);
            y = lpfs[idx].apply(u);
            EXPECT_LT(y, u);
        }
    }
}


/*
ResetTest:
- check, that the filter output value equals the input value after an reset (the filter should initialize with steady state on the current input value after an reset)
*/
TEST(LowPassFilterMpTest, ResetTest)
{
    static constexpr float sample_freq {1000};
    static constexpr float cutoff_freq {20};
    static constexpr uint8_t f_order {3};
    static constexpr uint8_t f_types[3] {LPF::FilterType::PTn, LPF::FilterType::Bessel, LPF::FilterType::Butterworth};
    
    static constexpr uint8_t sz = ARRAY_SIZE(f_types);

    const float u = 1.0F;
    float y;

    LPF lpfs[sz];
    for (uint8_t idx = 0; idx < sz; ++idx) {
        lpfs[idx].init(f_order, f_types[idx]);
        lpfs[idx].set_cutoff_frequency(sample_freq, cutoff_freq);
        lpfs[idx].reset();

        for (uint32_t i=0; i<100; ++i) {
            y = lpfs[idx].apply(u);
            
            EXPECT_NEAR(y, u, 1e-6F);
        }
    }
}


/*
ConfigChangeTest:
- check, that a change in sample_freq or cutoff_freq at runtime does not reset the filter
*/


/*
StepResponseTest:
- check for time constant?
- check for oszillation?
*/
TEST(LowPassFilterMpTest, StepResponseTest)
{
    static constexpr float sample_freq {1000};
    static constexpr float cutoff_freq {20};
    static constexpr uint8_t f_order {3};
    static constexpr uint8_t f_types[3] {LPF::FilterType::PTn, LPF::FilterType::Bessel, LPF::FilterType::Butterworth};
    
    static constexpr uint8_t  nt = ARRAY_SIZE(f_types);
    static constexpr uint32_t ns = (10.0F * 2.0F/(2.0F*M_PI*cutoff_freq)) * sample_freq;

    const float u = 1.0F;
    float y[nt][ns] {};

    LPF lpfs[nt];

    for (uint8_t t_idx = 0; t_idx < nt; ++t_idx) {
        lpfs[t_idx].init(f_order, f_types[t_idx]);
        lpfs[t_idx].set_cutoff_frequency(sample_freq, cutoff_freq);
        lpfs[t_idx].reset();
        lpfs[t_idx].apply(0.0F);

        for (uint32_t s_idx=0; s_idx<ns; ++s_idx) {
            y[t_idx][s_idx] = lpfs[t_idx].apply(u);
        }
    }

    // PTn:
    const float gain         = 1.0F/std::sqrt(2.0F);          // -3.01 dB
    const float gain_per_pt1 = std::pow(gain, 1.0F/f_order);  // gain per PT1 stage (gain must be splitted between _filter_order stages)
    const float tc_per_pt1   = std::sqrt(1.0F/std::pow(gain_per_pt1, 2) - 1.0F) / (cutoff_freq*2*M_PI);
    const float tc           = f_order * tc_per_pt1;

    uint32_t idx = std::round(tc * sample_freq);
    EXPECT_NEAR(0.63F*u, y[0][idx], 0.05F);
    // check also, all values smaller / greater; no overshoot; use google-mockup?


    // Bessel:


    // Butterworth:
}

/*
SineResponseTest:
- check for attenuation
*/

/*
ParametersTest
- check that calculated parameters are correct
*/

/*
//  test with a sine input
TEST(NotchFilterTest, SineTest)
{
    NotchFilter<float> filter;
    const float test_freq = 47;
    const float attenuation_dB = 43;
    const float rate_hz = 2000;
    const double dt = 1.0 / rate_hz;
    const uint32_t period_samples = rate_hz / test_freq;
    const uint32_t periods = 1000;
    const float test_amplitude = 0.7;
    const float expected_ratio = powf(10, (attenuation_dB/2)/10.0);
    double integral1_in = 0;
    double integral1_out = 0;
    double integral2_in = 0;
    double integral2_out = 0;
    filter.init(rate_hz, test_freq, test_freq*0.5, attenuation_dB);
    filter.reset();
    for (uint32_t i=0; i<periods * period_samples; i++) {
        const double t = i * dt;
        const double sample = sin(test_freq * t * 2 * M_PI) * test_amplitude;
        const float v = filter.apply(sample);
        if (i >= 2*period_samples) {
            integral1_in += sample * dt;
            integral2_in += fabsf(sample) * dt;
            integral1_out += v * dt;
            integral2_out += fabsf(v) * dt;
        }
    }

    // we expect both absolute integrals to be zero
    EXPECT_LE(fabsf(integral1_in), 0.01);
    EXPECT_LE(fabsf(integral1_out), 0.01);

    // we expect the output abs integral to be smaller than input
    // integral by the attenuation
    const float ratio1 = integral2_in / integral2_out;
    ::printf("ratio1=%f expected_ratio=%f\n", ratio1, expected_ratio);
    const float err_pct = 100 * fabsf(ratio1 - expected_ratio) / ratio1;
    EXPECT_LE(err_pct, 1);
}


//  test attentuation versus frequency
//  This is a way to get a graph of the attenuation and phase lag for a complex filter setup
TEST(NotchFilterTest, HarmonicNotchTest)
{
    const uint8_t num_test_freq = 150;
    const uint8_t harmonics = 15;
    const uint8_t num_harmonics = __builtin_popcount(harmonics);
    const float base_freq = 46;
    const float bandwidth = base_freq/2;
    const float attenuation_dB = 60;
    // number of identical filters chained together, simulating
    // usage of per-motor notch filtering
    const uint8_t chained_filters = 8;
    const uint16_t rate_hz = 2000;
    const uint32_t samples = 50000;
    const float test_amplitude = 0.7;
    const double dt = 1.0 / rate_hz;

    bool double_notch = true;
    HarmonicNotchFilter<float> filters[num_test_freq][chained_filters] {};
    struct {
        double in;
        double out;
        double last_in;
        double last_out;
        uint32_t last_crossing;
        uint32_t total_lag_samples;
        uint32_t lag_count;
        float get_lag_degrees(const float freq) const {
            const float lag_avg = total_lag_samples/float(lag_count);
            return (360.0 * lag_avg * freq) / rate_hz;
        }
    } integrals[num_test_freq] {};

    for (uint8_t i=0; i<num_test_freq; i++) {
        for (uint8_t c=0; c<chained_filters; c++) {
            auto &f = filters[i][c];
            f.allocate_filters(num_harmonics, harmonics, double_notch?2:1);
            f.init(rate_hz, base_freq, bandwidth, attenuation_dB);
        }
    }

    for (uint32_t s=0; s<samples; s++) {
        const double t = s * dt;

        for (uint8_t i=0; i<num_test_freq; i++) {
            const float freq = i+1;
            const double sample = sin(freq * t * 2 * M_PI) * test_amplitude;
            float v = sample;
            for (uint8_t c=0; c<chained_filters; c++) {
                auto &f = filters[i][c];
                v = f.apply(v);
            }
            if (s >= s/10) {
                integrals[i].in += fabsf(sample) * dt;
                integrals[i].out += fabsf(v) * dt;
            }
            if (sample >= 0 && integrals[i].last_in < 0) {
                integrals[i].last_crossing = s;
            }
            if (v >= 0 && integrals[i].last_out < 0 && integrals[i].last_crossing != 0) {
                integrals[i].total_lag_samples += s - integrals[i].last_crossing;
                integrals[i].lag_count++;
            }
            integrals[i].last_in = sample;
            integrals[i].last_out = v;
        }
    }
    const char *csv_file = "harmonicnotch_test.csv";
    FILE *f = fopen(csv_file, "w");
    fprintf(f, "Freq(Hz),Ratio,Lag(deg)\n");
    for (uint8_t i=0; i<num_test_freq; i++) {
        const float freq = i+1;
        const float lag_degrees = integrals[i].get_lag_degrees(freq);
        fprintf(f, "%.1f,%f,%f\n", freq, integrals[i].out/integrals[i].in, lag_degrees);
    }
    fclose(f);
    printf("Wrote %s\n", csv_file);
    ::printf("Lag at 1Hz %.2f degrees\n", integrals[0].get_lag_degrees(1));
    ::printf("Lag at 2Hz %.2f degrees\n", integrals[1].get_lag_degrees(2));
    ::printf("Lag at 5Hz %.2f degrees\n", integrals[4].get_lag_degrees(5));
    ::printf("Lag at 10Hz %.2f degrees\n", integrals[9].get_lag_degrees(10));
    EXPECT_NEAR(integrals[0].get_lag_degrees(1), 11.0, 0.5);
    EXPECT_NEAR(integrals[1].get_lag_degrees(2), 22.03, 0.5);
    EXPECT_NEAR(integrals[4].get_lag_degrees(5), 55.23, 0.5);
    EXPECT_NEAR(integrals[9].get_lag_degrees(10), 112.23, 0.5);
}
*/

AP_GTEST_MAIN()
