#pragma once

#include "arcscope/FilmTypes.h"

#include <sqlite3.h>

#include <string>

namespace arcscope {

class SQLiteStore {
  public:
    explicit SQLiteStore(const std::string& path);
    ~SQLiteStore();

    SQLiteStore(const SQLiteStore&) = delete;
    SQLiteStore& operator=(const SQLiteStore&) = delete;

    void initialize();
    void save_analysis(const FilmAnalysisResult& result);

  private:
    sqlite3* db_{nullptr};
};

}  // namespace arcscope

