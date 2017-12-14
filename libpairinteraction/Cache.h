/*
 * Copyright (c) 2016 Sebastian Weber, Henri Menke. All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef CACHE_H
#define CACHE_H

#include <mutex>
#include <unordered_map>

#include <boost/optional.hpp>


/** \brief Generic cache object
 *
 * This generic cache object strives to provide a thread-safe cache
 * based on a `std::unordered_map` protected from smashing by a mutex.
 *
 * There are only few public methods, which are save, restore, and
 * clear.  Iteration over the cache was left out intentionally as the
 * user should not search the cache manually for entries.
 */
template <typename Key, typename Element, typename Hash = std::hash<Key>>
class Cache
{
    typedef std::unordered_map<Key, Element, Hash> cache_t;

    cache_t cache;
    std::mutex cache_mutex;

public:
    /** \brief Save something in the cache
     *
     * To prevent race-conditions this function will throw if the
     * element is already in the cache.
     *
     * The usage is extremely straight-forward
     * \code
     * cache.save(key,element);
     * \endcode
     *
     * \param key Key
     * \param e Element
     * \throws std::runtime_error if the element is already in the cache
     */
    void save(Key const &key, Element const &e)
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        if (!cache.emplace(std::make_pair(key, e)).second) {
            throw std::runtime_error("Cache smashing detected!");
        }
    }

    /** \brief Restore something from the cache
     *
     * When restoring from the cache, the situation might occur that
     * an element is not yet available in the cache.  We do not want
     * to throw an exception in that case because that case might be
     * very frequent when building the cache.  Therefore restoring
     * returns an optional value, i.e. it might not void.  Hence the
     * access pattern is also a little different.
     *
     * \code
     * if (auto optional_element = cache.restore(key,element)) {
     *     // Cache hit!  Restore cache...
     *     element = optional_element.get();
     * } else {
     *     // Cache miss :-(
     * }
     * \endcode
     *
     * \param key Key
     * \returns Optional element
     */
    boost::optional<Element> restore(Key const &key)
    {
        std::lock_guard<std::mutex> lock(cache_mutex);
        auto cached_it = cache.find(key);
        if (cached_it != cache.end()) {
            return cached_it->second;
        }
        return boost::none;
    }

    /** \brief Clear the cache
     *
     * Delete all elements in the cache
     */
    void clear() { cache.clear(); }
};

#endif // CACHE_H
