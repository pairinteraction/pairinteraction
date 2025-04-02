// SPDX-FileCopyrightText: 2024 Pairinteraction Developers
// SPDX-License-Identifier: LGPL-3.0-or-later

#define DOCTEST_CONFIG_IMPLEMENT

#include "pairinteraction/tools/run_unit_tests.hpp"

#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/utils/paths.hpp"
#include "pairinteraction/utils/streamed.hpp"
#include "pairinteraction/version.hpp"

#include <cstdlib>
#include <doctest/doctest.h>
#include <filesystem>
#include <httplib.h>
#include <mutex>
#include <spdlog/spdlog.h>

// Create a reporter for doctest that logs to spdlog
namespace doctest {

// The code of the LoggingReporter is based on the ConsoleReporter from doctest,
// https://github.com/doctest/doctest/blob/ae7a13539fb71f270b87eb2e874fbac80bc8dda2/doctest/parts/doctest.cpp#L2868.
//
// SPDX-SnippetBegin
// SPDX-FileCopyrightText: (c) 2016-2025 Viktor Kirilov, Sebastian Weber
// SPDX-License-Identifier: MIT

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
                SPDLOG_ERROR(line);
            }
        } else {
            for (std::string line; std::getline(ss, line);) {
                SPDLOG_INFO(line);
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
            SPDLOG_INFO(ss.str());
        }

        if ((st.failure_flags & TestCaseFailureReason::Timeout) != 0) {
            std::stringstream ss;
            ss << Color::Red << "Test case exceeded time limit of " << std::setprecision(6)
               << std::fixed << tc->m_timeout << "!" << Color::None;
            SPDLOG_ERROR(ss.str());
        }

        if ((st.failure_flags & TestCaseFailureReason::ShouldHaveFailedButDidnt) != 0) {
            std::stringstream ss;
            ss << Color::Red << "Should have failed but didn't! Marking it as failed!"
               << Color::None;
            SPDLOG_ERROR(ss.str());
        } else if ((st.failure_flags & TestCaseFailureReason::ShouldHaveFailedAndDid) != 0) {
            std::stringstream ss;
            ss << Color::Yellow << "Failed as expected so marking it as not failed" << Color::None;
            SPDLOG_WARN(ss.str());
        } else if ((st.failure_flags & TestCaseFailureReason::CouldHaveFailedAndDid) != 0) {
            std::stringstream ss;
            ss << Color::Yellow << "Allowed to fail so marking it as not failed" << Color::None;
            SPDLOG_WARN(ss.str());
        } else if ((st.failure_flags & TestCaseFailureReason::DidntFailExactlyNumTimes) != 0) {
            std::stringstream ss;
            ss << Color::Red << "Didn't fail exactly " << tc->m_expected_failures
               << " times so marking it as failed!" << Color::None;
            SPDLOG_ERROR(ss.str());
        } else if ((st.failure_flags & TestCaseFailureReason::FailedExactlyNumTimes) != 0) {
            std::stringstream ss;
            ss << Color::Yellow << "Failed exactly " << tc->m_expected_failures
               << " times as expected so marking it as not failed!" << Color::None;
            SPDLOG_WARN(ss.str());
        }

        if ((st.failure_flags & TestCaseFailureReason::TooManyFailedAsserts) != 0) {
            std::stringstream ss;
            ss << Color::Red << "Aborting - too many failed asserts!" << Color::None;
            SPDLOG_ERROR(ss.str());
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
            SPDLOG_ERROR(loc + line);
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
                SPDLOG_ERROR(loc + line);
            }
        } else {
            for (std::string line; std::getline(ss, line);) {
                SPDLOG_INFO(loc + line);
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
            SPDLOG_INFO(loc + line);
        }

        log_contexts();
    }
};

// SPDX-SnippetEnd

REGISTER_REPORTER("logging", 1, doctest::LoggingReporter);
} // namespace doctest

constexpr std::string_view OS_NAME =
#if defined(_WIN32)
    "Windows";
#elif defined(__APPLE__)
    "macOS";
#elif defined(__linux__)
    "Linux";
#else
    "Unknown";
#endif

namespace pairinteraction {
int run_unit_tests(int argc, char **argv, bool download_missing, bool use_cache,
                   std::filesystem::path database_dir) {

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

    // Log the version and system information
    SPDLOG_INFO("Version of pairinteraction: {}.{}.{}", VERSION_MAJOR, VERSION_MINOR,
                VERSION_PATCH);
    SPDLOG_INFO("Operating system: {}", OS_NAME);

    // Create a global database instance and run the tests
    Database::get_global_instance(download_missing, use_cache, std::move(database_dir));
    int exitcode = ctx.run();

    std::filesystem::path logdir = paths::get_cache_directory() / "logs";
    SPDLOG_INFO("The log was stored to {}", logdir.string());

    if (exitcode != 0) {
        if (download_missing) {
            httplib::Client client("https://www.github.com");
            auto res = client.Head("/");
            if (!res) {
                SPDLOG_ERROR(
                    "Test failed. Please check your internet connection. An internet "
                    "connection is required to download databases of atomic states and matrix "
                    "elements if they are not available locally. The log was stored to {}",
                    logdir.string());
            } else {
                SPDLOG_ERROR(
                    "Tests failed. Consider creating an issue on "
                    "https://github.com/pairinteraction/pairinteraction/issues, attaching the "
                    "log. The log was stored to {}",
                    logdir.string());
            }
        } else {
            SPDLOG_ERROR(
                "Tests failed. Consider creating an issue on "
                "https://github.com/pairinteraction/pairinteraction/issues, attaching the "
                "log. If the tests failed because of unavailable states or "
                "matrix elements, consider downloading missing databases by calling "
                "the test function with 'download_missing = true'. The log was stored to {}",
                logdir.string());
        }
    }

    return exitcode;
};
} // namespace pairinteraction
