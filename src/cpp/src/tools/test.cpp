#define DOCTEST_CONFIG_IMPLEMENT

#include "pairinteraction/tools/test.hpp"

#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/utils/paths.hpp"
#include "pairinteraction/utils/streamed.hpp"

#include <cstdlib>
#include <doctest/doctest.h>
#include <filesystem>
#include <httplib.h>
#include <mutex>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

// Create a reporter for doctest that logs to spdlog
namespace doctest {

/*
 * The code of the LoggingReporter is based on the ConsoleReporter from doctest,
 * https://github.com/doctest/doctest/blob/ae7a13539fb71f270b87eb2e874fbac80bc8dda2/doctest/parts/doctest.cpp#L2868.
 * The following license applies to the code of the LoggingReporter:
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2016-2024 Viktor Kirilov, Sebastian Weber
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

// NOLINTBEGIN(cppcoreguidelines-macro-usage)
#define DOCTEST_LOCK_MUTEX(name)                                                                   \
    std::lock_guard<std::mutex> DOCTEST_ANONYMOUS(DOCTEST_ANON_LOCK_)(name);
// NOLINTEND(cppcoreguidelines-macro-usage)

struct LoggingReporter : public ConsoleReporter {
    LoggingReporter(const ContextOptions &co) : ConsoleReporter(co) {}

    void log_contexts() {}

    void logTestStart() {}

    void test_run_end(const TestRunStats &p) override {
        if (opt.minimal && p.numTestCasesFailed == 0) {
            return;
        }

        std::stringstream ss;
        ss << Color::Yellow
           << "==============================================================================="
           << Color::None << "\n";
        ss << std::dec;

        auto totwidth =
            int(std::ceil(log10(static_cast<double>(std::max(p.numTestCasesPassingFilters,
                                                             static_cast<unsigned>(p.numAsserts))) +
                                1)));
        auto passwidth =
            int(std::ceil(log10(static_cast<double>(std::max(
                                    p.numTestCasesPassingFilters - p.numTestCasesFailed,
                                    static_cast<unsigned>(p.numAsserts - p.numAssertsFailed))) +
                                1)));
        auto failwidth = int(
            std::ceil(log10(static_cast<double>(std::max(
                                p.numTestCasesFailed, static_cast<unsigned>(p.numAssertsFailed))) +
                            1)));
        const bool anythingFailed = p.numTestCasesFailed > 0 || p.numAssertsFailed > 0;
        ss << "test cases: " << std::setw(totwidth) << p.numTestCasesPassingFilters << " | "
           << ((p.numTestCasesPassingFilters == 0 || anythingFailed) ? Color::None : Color::Green)
           << std::setw(passwidth) << p.numTestCasesPassingFilters - p.numTestCasesFailed
           << " passed" << Color::None << " | "
           << (p.numTestCasesFailed > 0 ? Color::Red : Color::None) << std::setw(failwidth)
           << p.numTestCasesFailed << " failed" << Color::None << " |";
        if (!opt.no_skipped_summary) {
            const unsigned int numSkipped = p.numTestCases - p.numTestCasesPassingFilters;
            ss << " " << (numSkipped == 0 ? Color::None : Color::Yellow) << numSkipped << " skipped"
               << Color::None;
        }
        ss << "\n";
        ss << "assertions: " << std::setw(totwidth) << p.numAsserts << " | "
           << ((p.numAsserts == 0 || anythingFailed) ? Color::None : Color::Green)
           << std::setw(passwidth) << (p.numAsserts - p.numAssertsFailed) << " passed"
           << Color::None << " | " << (p.numAssertsFailed > 0 ? Color::Red : Color::None)
           << std::setw(failwidth) << p.numAssertsFailed << " failed" << Color::None << " |\n";
        ss << "Status: " << (p.numTestCasesFailed > 0 ? Color::Red : Color::Green)
           << ((p.numTestCasesFailed > 0) ? "FAILURE!" : "SUCCESS!") << Color::None << std::endl;

        if (p.numTestCasesFailed > 0) {
            for (std::string line; std::getline(ss, line);) {
                spdlog::error(line);
            }
        } else {
            for (std::string line; std::getline(ss, line);) {
                spdlog::info(line);
            }
        }
    }

    void test_case_end(const CurrentTestCaseStats &st) override {
        if (tc->m_no_output) {
            return;
        }

        if (opt.duration ||
            (st.failure_flags != 0 &&
             st.failure_flags != static_cast<int>(TestCaseFailureReason::AssertFailure))) {
            logTestStart();
        }

        if (opt.duration) {
            std::stringstream ss;
            ss << std::setprecision(6) << std::fixed << st.seconds << " s: " << tc->m_name;
            spdlog::info(ss.str());
        }

        if ((st.failure_flags & TestCaseFailureReason::Timeout) != 0) {
            std::stringstream ss;
            ss << Color::Red << "Test case exceeded time limit of " << std::setprecision(6)
               << std::fixed << tc->m_timeout << "!" << Color::None;
            spdlog::error(ss.str());
        }

        if ((st.failure_flags & TestCaseFailureReason::ShouldHaveFailedButDidnt) != 0) {
            std::stringstream ss;
            ss << Color::Red << "Should have failed but didn't! Marking it as failed!"
               << Color::None;
            spdlog::error(ss.str());
        } else if ((st.failure_flags & TestCaseFailureReason::ShouldHaveFailedAndDid) != 0) {
            std::stringstream ss;
            ss << Color::Yellow << "Failed as expected so marking it as not failed" << Color::None;
            spdlog::warn(ss.str());
        } else if ((st.failure_flags & TestCaseFailureReason::CouldHaveFailedAndDid) != 0) {
            std::stringstream ss;
            ss << Color::Yellow << "Allowed to fail so marking it as not failed" << Color::None;
            spdlog::warn(ss.str());
        } else if ((st.failure_flags & TestCaseFailureReason::DidntFailExactlyNumTimes) != 0) {
            std::stringstream ss;
            ss << Color::Red << "Didn't fail exactly " << tc->m_expected_failures
               << " times so marking it as failed!" << Color::None;
            spdlog::error(ss.str());
        } else if ((st.failure_flags & TestCaseFailureReason::FailedExactlyNumTimes) != 0) {
            std::stringstream ss;
            ss << Color::Yellow << "Failed exactly " << tc->m_expected_failures
               << " times as expected so marking it as not failed!" << Color::None;
            spdlog::warn(ss.str());
        }

        if ((st.failure_flags & TestCaseFailureReason::TooManyFailedAsserts) != 0) {
            std::stringstream ss;
            ss << Color::Red << "Aborting - too many failed asserts!" << Color::None;
            spdlog::error(ss.str());
        }
    }

    void test_case_exception(const TestCaseException &e) override {
        if (tc->m_no_output) {
            return;
        }

        DOCTEST_LOCK_MUTEX(mutex)

        logTestStart();

        std::stringstream ss;
        ss << "[" << skipPathFromFilename(tc->m_file.c_str()) << (opt.gnu_file_line ? ":" : "(")
           << (opt.no_line_numbers ? 0 : tc->m_line) << (opt.gnu_file_line ? "" : ")") << "] ";
        std::string loc = ss.str();
        ss.str("");
        ss << Color::Red << (e.is_crash ? "test case CRASHED: " : "test case THREW exception: ")
           << Color::None << e.error_string;
        for (std::string line; std::getline(ss, line);) {
            spdlog::error(loc + line);
        }
    }

    void log_assert(const AssertData &rb) override {
        if ((!rb.m_failed && !opt.success) || tc->m_no_output) {
            return;
        }

        DOCTEST_LOCK_MUTEX(mutex)

        logTestStart();

        std::stringstream ss;
        ss << "[" << skipPathFromFilename(rb.m_file) << (opt.gnu_file_line ? ":" : "(")
           << (opt.no_line_numbers ? 0 : rb.m_line) << (opt.gnu_file_line ? "" : ")") << "] ";
        std::string loc = ss.str();
        ss.str("");
        fulltext_log_assert_to_stream(ss, rb);
        if (rb.m_failed) {
            for (std::string line; std::getline(ss, line);) {
                spdlog::error(loc + line);
            }
        } else {
            for (std::string line; std::getline(ss, line);) {
                spdlog::info(loc + line);
            }
        }

        log_contexts();
    }

    void log_message(const MessageData &mb) override {
        if (tc->m_no_output) {
            return;
        }

        DOCTEST_LOCK_MUTEX(mutex)

        logTestStart();

        std::stringstream ss;
        ss << "[" << skipPathFromFilename(mb.m_file) << (opt.gnu_file_line ? ":" : "(")
           << (opt.no_line_numbers ? 0 : mb.m_line) << (opt.gnu_file_line ? "" : ")") << "] ";
        std::string loc = ss.str();
        ss.str("");
        ss << getSuccessOrFailColor(false, mb.m_severity)
           << getSuccessOrFailString((mb.m_severity & assertType::is_warn) != 0, mb.m_severity,
                                     "MESSAGE")
           << ": " << Color::None << mb.m_string;
        for (std::string line; std::getline(ss, line);) {
            spdlog::info(loc + line);
        }

        log_contexts();
    }
};

REGISTER_REPORTER("logging", 1, doctest::LoggingReporter);
} // namespace doctest

namespace pairinteraction {
int test(int argc, char **argv, bool download_missing, bool wigner_in_memory,
         std::filesystem::path database_dir) {

    // Configure a logger for the tests
    std::filesystem::path logdir = paths::get_pairinteraction_cache_directory() / "logs";

    if (!std::filesystem::exists(logdir)) {
        std::filesystem::create_directories(logdir);
    } else if (!std::filesystem::is_directory(logdir)) {
        throw std::runtime_error("Log path is not a directory.");
    }

    std::filesystem::path logfile = logdir / "test.log";

    auto console_sink =
        std::make_shared<spdlog::sinks::stdout_color_sink_mt>(spdlog::color_mode::always);
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(logfile.string(), true);
    auto doctest_logger =
        std::make_shared<spdlog::logger>(spdlog::logger("doctest", {console_sink, file_sink}));
    doctest_logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e %t] [%^%l%$] %v");

    const char *log_level = std::getenv("SPDLOG_LEVEL");
    if (log_level != nullptr) {
        spdlog::set_level(spdlog::level::from_str(log_level));
    }

    // Setup the tests
    doctest::Context ctx;
    ctx.setOption("abort-after", 5);
    ctx.setOption("no-run", 0);
    ctx.setOption("duration", true);
    ctx.setOption("no-path-filenames", true);
    ctx.applyCommandLine(argc, argv);
    ctx.setOption("no-colors", true);
    ctx.setOption("no-breaks", true);
    ctx.setOption("reporters", "logging");
    ctx.setOption("no-intro", true);

    // Set the doctest logger as the default logger
    auto original_logger = spdlog::default_logger();
    spdlog::set_default_logger(doctest_logger);

    // Create a global database instance and run the tests
    Database::get_global_instance(download_missing, wigner_in_memory, std::move(database_dir));
    int exitcode = ctx.run();

    spdlog::info("Log: {}", logfile.string());

    if (exitcode != 0) {
        if (download_missing) {
            httplib::Client client("https://www.github.com");
            auto res = client.Head("/");
            if (!res) {
                spdlog::error(
                    "Test failed. Please check your internet connection. An internet "
                    "connection is required to download databases of atomic states and matrix "
                    "elements if they are not available locally.",
                    logfile.string());
            } else {
                spdlog::error(
                    "Tests failed. Consider creating an issue on "
                    "https://github.com/pairinteraction/pairinteraction/issues, attaching the "
                    "log.",
                    logfile.string());
            }
        } else {
            spdlog::error(
                "Tests failed. Consider creating an issue on "
                "https://github.com/pairinteraction/pairinteraction/issues, attaching the "
                "log. If the tests failed because of unavailable states or "
                "matrix elements, consider downloading missing databases by calling "
                "the test function with 'download_missing = true'.",
                logfile.string());
        }
    }

    // Reset the default logger
    spdlog::set_default_logger(original_logger);

    return exitcode;
};
} // namespace pairinteraction
