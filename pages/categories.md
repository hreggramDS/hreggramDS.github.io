---
layout: page
title: "Posts by Category"
description: "Browse blog posts organized by category"
permalink: /blog/categories/
---

# Posts by Category

{% assign categories = site.categories | sort %}

{% for category in categories %}
<div class="category-section" id="{{ category[0] | slugify }}">

## {{ category[0] | capitalize }}

<span class="post-count">{{ category[1].size }} post{% if category[1].size != 1 %}s{% endif %}</span>

<ul class="category-posts">
{% for post in category[1] %}
  <li>
    <a href="{{ post.url | relative_url }}">{{ post.title }}</a>
    <span class="post-date">{{ post.date | date: "%B %d, %Y" }}</span>
  </li>
{% endfor %}
</ul>

</div>

---

{% endfor %}

{% if site.categories.size == 0 %}
<p><em>No categories yet. Categories will appear here as blog posts are published.</em></p>
{% endif %}
