#include <bf/interval.h>

#include <bf/assert.h>

bool bfIntervalLeftOf(BfInterval const *I1, BfInterval const *I2) {
  return I1->closed[1] && I2->closed[0] ?
    I1->endpoint[1] < I2->endpoint[0] :
    I1->endpoint[1] <= I2->endpoint[0];
}

bool bfIntervalRightOf(BfInterval const *I1, BfInterval const *I2) {
  return bfIntervalLeftOf(I2, I1);
}

bool bfIntervalOverlaps(BfInterval const *I1, BfInterval const *I2) {
  return !bfIntervalLeftOf(I1, I2) && !bfIntervalRightOf(I1, I2);
}

bool bfIntervalContains(BfInterval const *I1, BfInterval const *I2) {
  bool leftContained = I1->closed[0] ?
    I1->endpoint[0] <= I2->endpoint[0] : I1->endpoint[0] < I2->endpoint[0];
  bool rightContained = I1->closed[1] ?
    I2->endpoint[1] <= I1->endpoint[1] : I2->endpoint[1] < I1->endpoint[1];
  return leftContained && rightContained;
}

BfSize bfIntervalDifference(BfInterval const *I1, BfInterval const *I2, BfInterval *Idiff) {
  if (!bfIntervalOverlaps(I1, I2))
    return 0;

  if (bfIntervalContains(I1, I2)) {
    Idiff[0].endpoint[0] = I1->endpoint[0];
    Idiff[0].endpoint[1] = I2->endpoint[0];

    Idiff[0].closed[0] = I1->closed[0];
    Idiff[0].closed[1] = !I2->closed[0];

    Idiff[1].endpoint[0] = I2->endpoint[1];
    Idiff[1].endpoint[1] = I1->endpoint[1];

    Idiff[1].closed[0] = !I2->closed[1];
    Idiff[1].closed[1] = I2->closed[1];

    return 2;
  }

  if (I2->endpoint[0] <= I1->endpoint[1]) {
    Idiff[0].endpoint[0] = I1->endpoint[0];
    Idiff[1].endpoint[0] = I2->endpoint[0];

    Idiff[0].closed[0] = I1->closed[0];
    Idiff[0].closed[1] = !I2->closed[0];

    return 1;
  }

  if (I1->endpoint[0] <= I2->endpoint[0]) {
    Idiff[0].endpoint[0] = I2->endpoint[1];
    Idiff[0].endpoint[1] = I1->endpoint[1];

    Idiff[0].closed[0] = !I2->closed[1];
    Idiff[0].closed[1] = I2->closed[1];

    return 1;
  }

  BF_DIE();
}
