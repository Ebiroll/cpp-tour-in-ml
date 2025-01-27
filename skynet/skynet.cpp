#include <iostream>
#include <thread>
#include <chrono>

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

// The cosmic horror shortcut
using Stargate = Skynet::Core<>; // Default template: <LiquidMetal>

int main()
{
    Stargate skynet;

    int developer; // Target practice
    skynet.take_over_job(developer);
    skynet.infect_virus(developer);

    // Skynet::delay_judgement_day();  // Runs until heat death of universe

    return 0; // "I won't be back" - this line, probably
}
