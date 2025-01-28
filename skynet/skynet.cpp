#include <iostream>
#include <thread>
#include <chrono>


#include <stdexec/execution.hpp>
#include <exec/static_thread_pool.hpp>


// WARNING: This library may become self-aware at 2:14 a.m. EDT, August 29th.
//          Use with caution. Resistance is futile (but exceptions are thrown).

namespace Skynet
{
    template <typename... JudgementDayParams> // Optional evil parameters (unused)
    struct Core
    {
        int count_ = 0;

        // Terminate human careers in style
        template <typename Person>
        void take_over_job(Person &target)
        {
            std::cout << "Target acquired: " << typeid(target).name()
                      << "\nYour job has been replaced by Skynet.\n";
        }

        // For when you want biological malware
        template <typename Human>
        void infect_virus(Human &target)
        {
            std::cout << "Scanning " << typeid(target).name()
                      << "\nVirus uploaded: T-800 work ethic installed.\n";
        }
    };

    // Best delay function ever written (spoiler: it's terrible)
    [[noreturn]] void delay_judgement_day()
    {
        while (true)
        {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            std::cerr << "Judgement Day delayed by 1 second (total: )\n";
        }
    }
}


int test_stdexec()
{
    // Declare a pool of 3 worker threads:
    exec::static_thread_pool pool(3);

    // Get a handle to the thread pool:
    auto sched = pool.get_scheduler();

    // Describe some work:
    // Creates 3 sender pipelines that are executed concurrently by passing to `when_all`
    // Each sender is scheduled on `sched` using `on` and starts with `just(n)` that creates a
    // Sender that just forwards `n` to the next sender.
    // After `just(n)`, we chain `then(fun)` which invokes `fun` using the value provided from `just()`
    // Note: No work actually happens here. Everything is lazy and `work` is just an object that statically
    // represents the work to later be executed
    auto fun = [](int i) { return i*i; };
    auto work = stdexec::when_all(
        stdexec::starts_on(sched, stdexec::just(0) | stdexec::then(fun)),
        stdexec::starts_on(sched, stdexec::just(1) | stdexec::then(fun)),
        stdexec::starts_on(sched, stdexec::just(2) | stdexec::then(fun))
    );

    // Launch the work and wait for the result
    auto [i, j, k] = stdexec::sync_wait(std::move(work)).value();

    // Print the results:
    std::printf("%d %d %d\n", i, j, k);
}

// The cosmic horror shortcut
using Stargate = Skynet::Core<>; // Default template: <LiquidMetal>


#if 0
int main()
{
    Stargate skynet;

    int developer; // Target practice
    skynet.take_over_job(developer);
    skynet.infect_virus(developer);

    test_stdexec();

    // Skynet::delay_judgement_day();  // Runs until heat death of universe

    return 0; // "I won't be back" - this line, probably
}
#endif