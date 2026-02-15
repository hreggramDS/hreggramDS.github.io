source "https://rubygems.org"

# GitHub Pages gem - this includes Jekyll and all supported plugins
gem "github-pages", "~> 231", group: :jekyll_plugins

# Additional plugins (these are included in github-pages but explicitly listed for clarity)
group :jekyll_plugins do
  gem "jekyll-feed", "~> 0.12"
  gem "jekyll-seo-tag", "~> 2.8"
  gem "jekyll-sitemap", "~> 1.4"
  gem "jekyll-paginate", "~> 1.1"
end

# Required for Ruby 3.0+ (no longer a default gem)
gem "webrick", "~> 1.8"

# Silence CSV warning for Ruby 3.4+
gem "csv"

# Platform-specific gems for Windows and JRuby users
platforms :mingw, :x64_mingw, :mswin, :jruby do
  gem "tzinfo", ">= 1", "< 3"
  gem "tzinfo-data"
end

# Performance-booster for watching directories on Windows
gem "wdm", "~> 0.1", :platforms => [:mingw, :x64_mingw, :mswin]
