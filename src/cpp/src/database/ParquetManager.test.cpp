#include "pairinteraction/database/ParquetManager.hpp"

#include "pairinteraction/database/GitHubDownloader.hpp"

#include <doctest/doctest.h>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

namespace pairinteraction {
class MockDownloader : public GitHubDownloader {
public:
    explicit MockDownloader() : GitHubDownloader() {}

    std::future<GitHubDownloader::Result> download(const std::string &remote_url,
                                                   const std::string & /*if_modified_since*/ = "",
                                                   bool /*use_octet_stream*/ = false) override {
        GitHubDownloader::Result result;
        result.success = true;
        result.status_code = 200;
        result.rate_limit.remaining = 60;
        result.rate_limit.reset_time = 2147483647;

        if (remote_url == "/test/repo/path") {
            // This is the repo path request, return JSON with assets
            nlohmann::json assets = nlohmann::json::array();
            nlohmann::json asset;
            asset["name"] = "wigner_v3.parquet";
            asset["url"] = "https://api.github.com/test/path/wigner_v3.parquet";
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

TEST_CASE("ParquetManager basic functionality with mocked downloader") {
    MockDownloader downloader;
    auto test_dir = std::filesystem::temp_directory_path() / "pairinteraction_test_db";
    std::filesystem::create_directories(test_dir);
    std::ofstream(test_dir / "wigner_v1.parquet").close();
    std::ofstream(test_dir / "wigner_v2.parquet").close();

    SUBCASE("Check missing table") {
        std::vector<std::string> repo_paths;
        ParquetManager manager(test_dir, downloader, repo_paths);
        manager.scan_local();
        manager.scan_remote();

        CHECK_THROWS_WITH_AS(manager.get_path("missing_table"), "Table missing_table not found.",
                             std::runtime_error);
    }

    SUBCASE("Check version parsing") {
        std::vector<std::string> repo_paths;
        ParquetManager manager(test_dir, downloader, repo_paths);
        manager.scan_local();
        manager.scan_remote();

        std::string expected = (test_dir / "wigner_v2.parquet").string();
        CHECK(manager.get_path("wigner") == expected);
    }

    SUBCASE("Check update table") {
        std::vector<std::string> repo_paths = {"/test/repo/path"};
        std::ofstream(test_dir / "wigner_v1.parquet").close();
        std::ofstream(test_dir / "wigner_v2.parquet").close();
        ParquetManager manager(test_dir, downloader, repo_paths);
        manager.scan_local();
        manager.scan_remote();

        std::string expected = (test_dir / "wigner_v3.parquet").string();
        CHECK(manager.get_path("wigner") == expected);

        std::ifstream in(expected, std::ios::binary);
        std::stringstream buffer;
        buffer << in.rdbuf();
        CHECK(buffer.str() == "updated_file_content");
    }

    std::filesystem::remove_all(test_dir);
}

DOCTEST_TEST_CASE(
    "ParquetManager basic functionality with github downloader" *
    doctest::skip(true)) { // TODO skip if !Database::get_global_instance().download_missing()
    GitHubDownloader downloader;
    auto test_dir = std::filesystem::temp_directory_path() / "pairinteraction_test_db";

    std::vector<std::string> repo_paths = {"/repos/pairinteraction/database-sqdt/releases/latest",
                                           "/repos/pairinteraction/database-mqdt/releases/latest"};
    ParquetManager manager(test_dir, downloader, repo_paths);
    manager.scan_local();
    manager.scan_remote();

    CHECK(manager.get_path("Rb_states").rfind(test_dir / "Rb_states_v", 0) == 0);
    // TODO check only that the database overview has been downloaded correctly. Downloading a full
    // table takes too long.

    std::filesystem::remove_all(test_dir);
}

} // namespace pairinteraction
