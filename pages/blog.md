---
layout: page
title: Blog
description: Technical posts, tutorials, and research updates
permalink: /blog/
---

# Blog

Technical writing on data science, machine learning, bioinformatics, and scientific computing.

## Recent Posts

{% assign posts = site.posts | sort: "date" | reverse %}

{% if posts.size > 0 %}
<div class="blog-posts">
  {% for post in posts limit:10 %}
  <article class="blog-post-preview">
    <h2><a href="{{ post.url | relative_url }}">{{ post.title }}</a></h2>
    <div class="post-meta">
      <time datetime="{{ post.date | date_to_xmlschema }}">
        {{ post.date | date: "%B %d, %Y" }}
      </time>
      {% if post.categories %}
      <span class="post-categories">
        in {% for category in post.categories %}
          <a href="{{ '/blog/categories/' | append: category | relative_url }}">{{ category }}</a>{% unless forloop.last %}, {% endunless %}
        {% endfor %}
      </span>
      {% endif %}
    </div>
    
    {% if post.excerpt %}
    <p class="post-excerpt">{{ post.excerpt | strip_html | truncatewords: 50 }}</p>
    {% endif %}
    
    {% if post.tags %}
    <div class="post-tags">
      {% for tag in post.tags limit:5 %}
        <span class="tag">{{ tag }}</span>
      {% endfor %}
    </div>
    {% endif %}
    
    <a href="{{ post.url | relative_url }}" class="read-more">Read more →</a>
  </article>
  <hr>
  {% endfor %}
</div>

{% if paginator.total_pages > 1 %}
<nav class="pagination">
  {% if paginator.previous_page %}
    <a href="{{ paginator.previous_page_path | relative_url }}" class="prev">← Newer Posts</a>
  {% endif %}
  
  <span class="page-number">Page {{ paginator.page }} of {{ paginator.total_pages }}</span>
  
  {% if paginator.next_page %}
    <a href="{{ paginator.next_page_path | relative_url }}" class="next">Older Posts →</a>
  {% endif %}
</nav>
{% endif %}

{% else %}
<p><em>Blog posts coming soon! Check back for technical tutorials, project deep-dives, and research updates.</em></p>

## Planned Topics

- **Single-cell RNA-seq Analysis**: Complete tutorial series
- **Topological Data Analysis**: Practical applications
- **Machine Learning in Production**: Best practices and lessons learned
- **From Physics to Data Science**: Career transition insights
- **Python for Scientific Computing**: Advanced techniques
{% endif %}

---

## Categories

{% if site.categories.size > 0 %}
<div class="category-list">
  {% for category in site.categories %}
    <a href="{{ '/blog/categories/' | append: category[0] | relative_url }}" class="category-badge">
      {{ category[0] }} ({{ category[1].size }})
    </a>
  {% endfor %}
</div>
{% endif %}

## Subscribe

*Optional: Add RSS feed link or newsletter signup*

Stay updated with new posts via [RSS](/feed.xml).

---

## Writing Goals

This blog serves as a platform for:
- **Technical Documentation**: Sharing methodologies and best practices
- **Tutorials**: Step-by-step guides for complex analyses
- **Project Narratives**: Stories behind the data and insights
- **Knowledge Sharing**: Lessons learned from real-world applications
