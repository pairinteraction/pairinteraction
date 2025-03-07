#include "pairinteraction/database/ParquetManager.hpp"

#include "pairinteraction/database/Database.hpp"
#include "pairinteraction/database/GitHubDownloader.hpp"

#include <doctest/doctest.h>
#include <duckdb.hpp>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

namespace pairinteraction {
class MockDownloader : public GitHubDownloader {
public:
    std::future<GitHubDownloader::Result>
    download(const std::string &remote_url, const std::string & /*if_modified_since*/ = "",
             bool /*use_octet_stream*/ = false) const override {
        GitHubDownloader::Result result;
        result.status_code = 200;
        result.rate_limit.remaining = 60;
        result.rate_limit.reset_time = 2147483647;

        if (remote_url == "/test/repo/path") {
            // This is the repo path request, return JSON with assets
            nlohmann::json assets = nlohmann::json::array();
            nlohmann::json asset;
            asset["name"] = "wigner_v1.2.parquet";
            asset["url"] = "https://api.github.com/test/path/wigner_v1.2.parquet";
            assets.push_back(asset);
            nlohmann::json response;
            response["assets"] = assets;
            result.body = response.dump();
        } else if (remote_url == "/rate_limit") {
            // This is the rate limit request
            result.body = "";
        } else {
            // This is the file download request
            result.body = "updated_file_content";
        }

        return std::async(std::launch::deferred, [result]() { return result; });
    }
};

TEST_CASE("ParquetManager functionality with mocked downloader") {
    MockDownloader downloader;
    auto test_dir = std::filesystem::temp_directory_path() / "pairinteraction_test_db";
    std::filesystem::create_directories(test_dir / "tables");
    std::ofstream(test_dir / "tables" / "wigner_v1.0.parquet").close();
    std::ofstream(test_dir / "tables" / "wigner_v1.1.parquet").close();
    duckdb::DuckDB db(nullptr);
    duckdb::Connection con(db);

    SUBCASE("Check missing table") {
        std::vector<std::string> repo_paths;
        ParquetManager manager(test_dir, downloader, repo_paths, con, false);
        manager.scan_local();
        manager.scan_remote();

        CHECK_THROWS_WITH_AS(manager.get_path("", "missing_table"),
                             "Table missing_table not found.", std::runtime_error);
    }

    SUBCASE("Check version parsing") {
        std::vector<std::string> repo_paths;
        ParquetManager manager(test_dir, downloader, repo_paths, con, false);
        manager.scan_local();
        manager.scan_remote();

        std::string expected = (test_dir / "tables" / "wigner_v1.1.parquet").string();
        CHECK(manager.get_path("", "wigner") == expected);
    }

    SUBCASE("Check update table") {
        std::vector<std::string> repo_paths = {"/test/repo/path"};
        ParquetManager manager(test_dir, downloader, repo_paths, con, false);
        manager.scan_local();
        manager.scan_remote();

        std::string expected = (test_dir / "tables" / "wigner_v1.2.parquet").string();
        CHECK(manager.get_path("", "wigner") == expected);

        std::ifstream in(expected, std::ios::binary);
        std::stringstream buffer;
        buffer << in.rdbuf();
        CHECK(buffer.str() == "updated_file_content");
    }

    std::filesystem::remove_all(test_dir);
}

DOCTEST_TEST_CASE("ParquetManager functionality with github downloader" * doctest::skip(true)) {
    // TODO This test is currently skipped because not all tables are yet available

    if (!Database::get_global_instance().get_download_missing()) {
        DOCTEST_MESSAGE("Skipping test because download_missing is false.");
        return;
    }
    GitHubDownloader downloader;
    duckdb::DuckDB db(nullptr);
    duckdb::Connection con(db);

    std::vector<std::string> repo_paths = {"/repos/pairinteraction/database-sqdt/releases/latest",
                                           "/repos/pairinteraction/database-mqdt/releases/latest"};
    ParquetManager manager(Database::get_global_instance().get_database_dir(), downloader,
                           repo_paths, con, Database::get_global_instance().get_use_cache());
    manager.scan_local();
    manager.scan_remote();

    std::string info = manager.get_versions_info();

    // Check that all species are present
    std::vector<std::string> should_contain = {
        "Cs",        "K",          "Li",           "Na",
        "Rb",        "Sr87_mqdt",  "Sr88_singlet", "Sr88_triplet",
        "Sr88_mqdt", "Yb171_mqdt", "Yb173_mqdt",   "Yb174_mqdt",
        "wigner"};
    for (const auto &substr : should_contain) {
        DOCTEST_CHECK(info.find(substr) != std::string::npos);
    }
}

} // namespace pairinteraction
