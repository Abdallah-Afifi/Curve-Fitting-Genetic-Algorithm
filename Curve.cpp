#include <bits/stdc++.h>
#define endl '\n'
#define ll long long
using namespace std;

const int popSize = 40, numIterations = 100;
const double crossProb = 0.5, mutProb = 0.05;
int countGen = 0;

class Chromosome {
private:
    ll degree = 0;
    long double fitness{};
    long double totalErr{};
    vector<pair<long double, long double>> allItems;
    int itemNum{};

public:
    vector<long double> genes;
    Chromosome() {
        itemNum = ++countGen;
    }
    virtual ~Chromosome() = default;
    void setItems(ll d) {
        this->degree = d;
        genes.resize(this->degree);
    }
    void calcFitness(const vector<pair<long double, long double>>& items) {
        allItems = items;
        long double sum, finalErr = 0.0;
        for (const auto& item : items) {
            sum = genes[0];
            for (int j = 1; j < genes.size(); j++)
                sum += (genes[j] * pow(item.first, j));
            sum -= item.second;
            sum *= sum;
            finalErr += sum;
        }
        this->totalErr = (finalErr / items.size());
        this->fitness = (1 / this->totalErr);
    }
    [[nodiscard]] long double getFitness() const {
        return fitness;
    }
    [[nodiscard]] long double getTotalErr() const {
        return totalErr;
    }
    [[nodiscard]] vector<pair<long double, long double>> getAllItems() const {
        return allItems;
    }
};

bool compareValues(const Chromosome& c1, const Chromosome& c2) {
    return c1.getFitness() > c2.getFitness();
}

class GeneticAlgorithm {
private:
    vector<Chromosome>* pop;
    vector<Chromosome> offspring;
    int numOffspring;

    Chromosome tournament() {
        auto& population = *pop;
        Chromosome bestChromosome, chromosome;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> distribution(0, (int)population.size() - 3);
        bestChromosome = population[distribution(gen)];

        for (int i = 1; i < this->numOffspring; i++) {
            chromosome = population[distribution(gen)];
            chromosome.calcFitness(population[0].getAllItems());
            bestChromosome.calcFitness(population[0].getAllItems());
            if (chromosome.getFitness() < bestChromosome.getFitness())
                bestChromosome = chromosome;
        }
        return bestChromosome;
    }

    void elitism() {
        auto& population = *pop;
        sort(population.begin(), population.end(), compareValues);
        for (int i = 0; i < offspring.size(); i++)
            population[i] = offspring[i];
    }

public:
    explicit GeneticAlgorithm(vector<Chromosome>& population) {
        this->pop = &population;
        this->numOffspring = (int)this->pop->size() / 2;
        offspring.resize(this->numOffspring);
    }

    void initializePopulation() {
        for (auto& individual : *pop) {
            for (auto&& gene : individual.genes) {
                gene = ((float(rand()) / float(RAND_MAX)) * (10 + 10)) - 10;
            }
        }
    }

    void selection() {
        auto& population = *pop;
        sort(population.begin(), population.end(), compareValues);
        for (int i = 0; i < numOffspring; i++)
            offspring[i] = this->tournament();
    }

    void crossOver() {
        for (int i = 0; i < offspring.size(); i++) {
            for (int j = i + 1; j < offspring.size(); j++) {
                double randNum = ((double)rand() / (RAND_MAX));
                if (randNum < crossProb) {
                    int randIdx = rand() % (offspring[i].genes.size() - 2);
                    for (int k = randIdx; k < offspring[i].genes.size(); k++)
                        swap(offspring[i].genes[k], offspring[j].genes[k]);
                }
            }
        }
    }

    void mutation(ll t) {
        for (auto& individual : offspring) {
            long double maxVal = -1e9, minVal = 1e9;
            for (auto&& gene : individual.genes)
                if (gene > maxVal) maxVal = gene;

            for (auto&& gene : individual.genes)
                if (gene < minVal) minVal = gene;

            for (auto&& gene : individual.genes) {
                double randNum = ((double)rand() / (RAND_MAX));
                if (randNum < mutProb) {
                    long double dLower = gene - minVal, dUpper = maxVal - gene, y = 0, newVal;
                    double randNum1 = ((double)rand() / (RAND_MAX));
                    y = randNum1 <= 0.5 ? dLower : dUpper;
                    newVal = y * (1 - pow(randNum1, pow(1 - t / numIterations, 2)));
                    gene = (y == dLower) ? gene - newVal : gene + newVal;
                }
            }
        }
    }

    void replacement() {
        elitism();
    }
};

int main() {
    srand(time(nullptr));
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    ll numTestCases, numItems, degree, iterSize, testNum = 1;
    cin >> numTestCases;
    while (numTestCases--) {
        iterSize = numIterations;
        cin >> numItems >> degree;
        vector<pair<long double, long double>> items(numItems);
        for (auto i = 0; i < numItems; i++)
            cin >> items[i].first >> items[i].second;
        vector<Chromosome> population(popSize);
        for (int i = 0; i < popSize; i++)
            population[i].setItems(degree + 1);
        GeneticAlgorithm geneticAlgorithm(population);

        // First Step (Initialize Pool Of Solutions)
        geneticAlgorithm.initializePopulation();

        while (iterSize--) {
            // Second Step (Individual Evaluation)
            for (int i = 0; i < popSize; i++)
                population[i].calcFitness(items);

            // Third Step (Selection)
            geneticAlgorithm.selection();

            // Fourth Step (Crossover)
            geneticAlgorithm.crossOver();

            // Fifth Step (Mutation)
            geneticAlgorithm.mutation(testNum);

            // Sixth Step (Reproduction)
            geneticAlgorithm.replacement();
        }
        sort(population.begin(), population.end(), compareValues);

        // Printing
        cout << "Testcase " << testNum++ << ":\nThe Coefficients: [ ";
        for (auto&& gene : population[population.size() - 1].genes) {
            cout << fixed << setprecision(2);
            gene != population[population.size() - 1].genes[population[population.size() - 1].genes.size() - 1] ?
                cout << gene << " , " :
                cout << gene << " ";
        }
        cout << "]" << endl << "Mean Square Error: " << fixed << setprecision(20)
             << population[population.size() - 1].getFitness() << endl << endl;
    }
    return 0;
}
